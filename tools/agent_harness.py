#!/usr/bin/env python3
"""
Multi-agent orchestrator: Copilot CLI (ACP) ↔ Kiro CLI (ACP).

Architecture
~~~~~~~~~~~~
  tools/agent_harness.py sits alongside the ACOLITE codebase as a
  development-time orchestration utility.  It is NOT part of the
  atmospheric-correction pipeline (acolite/ or src/).

  Repository layout:
    acolite/          ← Python ACOLITE library (48+ sensors)
    src/              ← Rust port (acolite-rs)
    tests/            ← Rust E2E + regression tests
    examples/         ← Rust processing examples
    config/           ← Settings / LUT metadata
    tools/            ← **Development tooling (this file)**
      agent_harness.py

Multi-Agent Diagram  (see RUST_PORT_ROADMAP.md for rendered Mermaid)
~~~~~~~~~~~~~~~~~~~
  User --task --> [Orchestrator]
      |                                |
      v                                v
  [Kiro CLI]  <--- ACP/stdio --->  [Copilot CLI]
  (Executor)                       (Proposer)
  kiro-cli acp                     copilot --acp --stdio
      |                                |
      +-----> Proposal Loop <----------+
                   |
                   v
            [Human Console]
            (approve/reject)

  Data flow:
    1. User --task --> Orchestrator --> Kiro (ACP prompt)
    2. Kiro streams output via session/update
    3. Kiro output --> Copilot (ACP prompt for review)
    4. Copilot proposes ACTION: lines
    5. Human approves --> approved actions --> Kiro
    6. Repeat 2-5 until DONE or max cycles

Protocol Notes
~~~~~~~~~~~~~~
Both CLIs support the Agent Client Protocol (ACP) over NDJSON stdio:
  - Copilot CLI: copilot --acp --stdio
    Ref: https://docs.github.com/en/copilot/reference/acp-server
  - Kiro CLI:    kiro-cli acp [--agent NAME]
    Ref: https://kiro.dev/docs/cli/acp/

ACP uses JSON-RPC 2.0 over stdin/stdout with these core methods:
  initialize → session/new → session/prompt → (session/update notifications)

Key differences between the two agents:
  - protocolVersion: integer 1 for both Copilot and Kiro
  - session/prompt field: Copilot "prompt" vs Kiro "content"
  - Kiro extensions: _kiro.dev/commands/*, _kiro.dev/mcp/*

Usage
~~~~~
    python tools/agent_harness.py --task "Port gas_transmittance to Rust"

    python tools/agent_harness.py \\
        --task "Fix failing tests" \\
        --kiro-agent rust-developer \\
        --auto-approve
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(name)-12s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("harness")


# ---------------------------------------------------------------------------
# ACP client — NDJSON over stdio, works with both Copilot and Kiro CLIs
# Ref: https://agentclientprotocol.com/protocol/overview
# ---------------------------------------------------------------------------

class ACPClient:
    """
    Communicates with an ACP agent over NDJSON (newline-delimited JSON-RPC 2.0)
    via stdio. Works with both Copilot CLI and Kiro CLI.

    Copilot: copilot --acp --stdio
    Kiro:    kiro-cli acp [--agent NAME]

    Both implement the same core ACP methods: initialize, session/new,
    session/prompt, session/cancel. Responses and session/update notifications
    arrive as NDJSON lines on stdout.

    agent_type: "copilot" or "kiro" — both now send protocolVersion as integer 1.
    """

    def __init__(self, process: subprocess.Popen, name: str = "agent",
                 agent_type: str = "kiro"):
        self.process = process
        self.name = name
        self.agent_type = agent_type  # "copilot" | "kiro"
        self._id_counter = 0
        self._lock = threading.Lock()
        self._pending: dict[str, threading.Event] = {}
        self._results: dict[str, dict] = {}
        self._session_id: Optional[str] = None
        self._text_lock = threading.Lock()   # guards _streamed_text
        self._streamed_text: list[str] = []
        self._tool_calls: list[dict] = []
        self._reader_thread = threading.Thread(
            target=self._read_loop, daemon=True
        )
        self._reader_thread.start()

    def _next_id(self) -> str:
        with self._lock:
            self._id_counter += 1
            return str(self._id_counter)

    def _send(self, method: str, params: dict) -> dict:
        """Send a JSON-RPC 2.0 request and wait for the response."""
        msg_id = self._next_id()
        payload = {
            "jsonrpc": "2.0",
            "id": msg_id,
            "method": method,
            "params": params,
        }
        event = threading.Event()
        with self._lock:
            self._pending[msg_id] = event

        line = json.dumps(payload) + "\n"
        try:
            self.process.stdin.write(line.encode())
            self.process.stdin.flush()
        except (BrokenPipeError, OSError) as e:
            log.error("[%s] Write failed: %s", self.name, e)
            return {"error": str(e)}

        if not event.wait(timeout=300):
            log.warning("[%s] Timeout waiting for response to %s", self.name, method)
            return {"error": "timeout"}

        return self._results.pop(msg_id, {"error": "no result"})

    def _notify(self, method: str, params: dict):
        """Send a JSON-RPC 2.0 notification (no response expected)."""
        payload = {
            "jsonrpc": "2.0",
            "method": method,
            "params": params,
        }
        line = json.dumps(payload) + "\n"
        try:
            self.process.stdin.write(line.encode())
            self.process.stdin.flush()
        except (BrokenPipeError, OSError) as e:
            log.error("[%s] Notify write failed: %s", self.name, e)

    def _read_loop(self):
        """Read NDJSON responses from the agent's stdout."""
        for raw_line in self.process.stdout:
            line = raw_line.decode().strip() if isinstance(raw_line, bytes) else raw_line.strip()
            if not line:
                continue
            try:
                msg = json.loads(line)
            except json.JSONDecodeError:
                log.debug("[%s] Non-JSON line: %s", self.name, line[:100])
                continue

            # JSON-RPC response (has "id") — matches a pending request
            msg_id = msg.get("id")
            if msg_id is not None:
                msg_id = str(msg_id)
                with self._lock:
                    if msg_id in self._pending:
                        self._results[msg_id] = msg.get("result", msg)
                        self._pending[msg_id].set()
                continue

            # Notification from agent (no "id" field)
            method = msg.get("method", "")
            params = msg.get("params", {})
            self._handle_notification(method, params)

    def _handle_notification(self, method: str, params: dict):
        """Handle ACP session/update notifications from the agent."""
        if method == "session/update":
            update = params.get("update", {})
            update_type = update.get("sessionUpdate", "")

            if update_type == "agent_message_chunk":
                content = update.get("content", {})
                if content.get("type") == "text":
                    text = content.get("text", "")
                    if text:
                        with self._text_lock:
                            self._streamed_text.append(text)
                        print(text, end="", flush=True)

            elif update_type == "tool_call":
                tool_info = {
                    "name": update.get("name", ""),
                    "status": update.get("status", ""),
                    "params": update.get("parameters", {}),
                }
                self._tool_calls.append(tool_info)
                log.debug("[%s] Tool call: %s (%s)",
                          self.name, tool_info["name"], tool_info["status"])

            elif update_type == "turn_end":
                log.debug("[%s] Turn ended", self.name)

        elif method == "session/request_permission":
            # Grant all requested tool permissions so agents can work unblocked.
            # ACP spec: respond with session/grant_permission notification
            # carrying the same permissionId(s) from the request.
            permissions = params.get("permissions", [])
            permission_ids = [p.get("permissionId") for p in permissions
                              if isinstance(p, dict) and p.get("permissionId")]
            log.info("[%s] Permission request — granting: %s",
                     self.name, permission_ids or permissions)
            self._notify("session/grant_permission", {
                "sessionId": params.get("sessionId") or self._session_id,
                "permissionIds": permission_ids,
                "granted": True,
            })

        # Kiro custom extensions — log but don't act
        elif method.startswith("_kiro.dev/"):
            log.debug("[%s] Kiro extension: %s", self.name, method)

    # -- ACP lifecycle methods ---------------------------------------------

    def initialize(self, client_name: str = "acolite-harness") -> dict:
        """ACP initialize — negotiate protocol version and capabilities.

        Both Copilot and Kiro use protocolVersion as an integer (1).
        """
        protocol_version: Any = 1
        return self._send("initialize", {
            "protocolVersion": protocol_version,
            "clientCapabilities": {
                "fs": {"readTextFile": True, "writeTextFile": True},
                "terminal": True,
            },
            "clientInfo": {"name": client_name, "version": "1.0.0"},
        })

    def new_session(self, cwd: str,
                    mcp_servers: Optional[list] = None) -> Optional[str]:
        """ACP session/new — create a coding session."""
        result = self._send("session/new", {
            "cwd": cwd,
            "mcpServers": mcp_servers or [],
        })
        self._session_id = result.get("sessionId")
        return self._session_id

    def prompt(self, text: str, session_id: Optional[str] = None) -> dict:
        """ACP session/prompt — send a user message.

        Both Copilot and Kiro accept content as an array of content parts.
        Copilot uses "prompt" key, Kiro uses "content" key.
        We send both for compatibility.
        """
        sid = session_id or self._session_id
        if not sid:
            log.error("[%s] No session ID — call new_session() first", self.name)
            return {"error": "no session"}
        content_parts = [{"type": "text", "text": text}]
        return self._send("session/prompt", {
            "sessionId": sid,
            "prompt": content_parts,    # Copilot expects "prompt"
            "content": content_parts,   # Kiro expects "content"
        })

    def cancel(self, session_id: Optional[str] = None):
        """ACP session/cancel — interrupt current processing."""
        sid = session_id or self._session_id
        if sid:
            self._notify("session/cancel", {"sessionId": sid})

    def drain_streamed_text(self) -> str:
        """Return and clear all text accumulated from session/update notifications."""
        with self._text_lock:
            text = "".join(self._streamed_text)
            self._streamed_text.clear()
        return text

    def drain_tool_calls(self) -> list[dict]:
        """Return and clear tool call history."""
        calls = list(self._tool_calls)
        self._tool_calls.clear()
        return calls

    def alive(self) -> bool:
        """Check if the underlying process is still running."""
        return self.process.poll() is None

    def shutdown(self):
        """Clean shutdown of the agent process."""
        try:
            if self.process.stdin:
                self.process.stdin.close()
            self.process.terminate()
            self.process.wait(timeout=5)
        except Exception:
            self.process.kill()


# ---------------------------------------------------------------------------
# Agent launchers
# ---------------------------------------------------------------------------

def launch_copilot_acp(workdir: str) -> subprocess.Popen:
    """Start Copilot CLI as an ACP agent: copilot --acp --stdio"""
    executable = os.environ.get("COPILOT_CLI_PATH", "copilot")
    cmd = [executable, "--acp", "--stdio"]
    log.info("Launching %s", " ".join(cmd))
    return subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=workdir,
    )


def launch_kiro_acp(workdir: str,
                    agent: Optional[str] = None) -> subprocess.Popen:
    """Start Kiro CLI as an ACP agent: kiro-cli acp [--agent NAME]

    Ref: https://kiro.dev/docs/cli/acp/
    """
    executable = os.environ.get("KIRO_CLI_PATH", "kiro-cli")
    cmd = [executable, "acp"]
    if agent:
        cmd.extend(["--agent", agent])
    log.info("Launching %s", " ".join(cmd))
    return subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=workdir,
    )


def init_acp_agent(proc: subprocess.Popen, name: str,
                   workdir: str, agent_type: str = "kiro") -> ACPClient:
    """Create an ACPClient, initialize connection, and create a session."""
    client = ACPClient(proc, name=name, agent_type=agent_type)

    init_resp = client.initialize(client_name="acolite-harness")
    if "error" in init_resp:
        client.shutdown()
        raise RuntimeError(f"{name} ACP initialize failed: {init_resp['error']}")

    agent_info = init_resp.get("agentInfo", {})
    log.info("[%s] Initialized: %s v%s",
             name, agent_info.get("name", "?"), agent_info.get("version", "?"))

    session_id = client.new_session(workdir)
    if not session_id:
        client.shutdown()
        raise RuntimeError(f"{name} session/new failed — no sessionId returned")

    log.info("[%s] Session: %s", name, session_id)
    return client


# ---------------------------------------------------------------------------
# Orchestrator — the routing loop
# ---------------------------------------------------------------------------

class Orchestrator:
    """
    Routes between Kiro (executor, ACP) and Copilot (proposer, ACP).
    Both agents communicate via the same ACP protocol over stdio.
    Kiro executes tasks, Copilot reviews progress and proposes next steps.
    """

    def __init__(self, kiro: ACPClient, copilot: ACPClient,
                 workdir: str, auto_approve: bool = False,
                 max_cycles: int = 10):
        self.kiro = kiro
        self.copilot = copilot
        self.workdir = workdir
        self.auto_approve = auto_approve
        self.max_cycles = max_cycles
        self.event_log: list[dict] = []
        self._kiro_output: list[str] = []

    def run(self, task_description: str) -> dict:
        """Main orchestration loop."""
        task_id = str(uuid.uuid4())[:12]

        # 1. Send initial task to Kiro via ACP session/prompt
        log.info("Sending task to Kiro: %s", task_description[:80])
        self._log_event("task_submitted", task_id, "kiro")
        self.kiro.drain_streamed_text()  # clear any prior text

        kiro_resp = self.kiro.prompt(task_description)
        kiro_text = self.kiro.drain_streamed_text()
        self._kiro_output.append(kiro_text)

        stop_reason = kiro_resp.get("stopReason", "")
        self._log_event("kiro_response", task_id, "kiro",
                        stop_reason=stop_reason, chars=len(kiro_text))
        print(flush=True)  # newline after streamed output

        if "error" in kiro_resp:
            log.error("Kiro prompt failed: %s", kiro_resp["error"])
            return {"status": "failed", "output": kiro_text,
                    "events": len(self.event_log)}

        # 2. Iterative proposal cycles
        self._proposal_loop(task_id, task_description)

        return {
            "status": "completed",
            "output": "\n---\n".join(self._kiro_output[-3:]),
            "events": len(self.event_log),
        }

    def _proposal_loop(self, task_id: str, original_task: str):
        """Copilot proposes, human approves, Kiro executes."""
        for cycle in range(1, self.max_cycles + 1):
            # Build context from Kiro's recent output
            context = "\n---\n".join(self._kiro_output[-5:])
            if not context.strip():
                break

            # Ask Copilot for proposals via ACP
            log.info("Proposal cycle %d", cycle)
            print(f"\n  [copilot] Analyzing progress (cycle {cycle})...",
                  flush=True)
            self.copilot.drain_streamed_text()  # clear before prompt

            self.copilot.prompt(
                f"An executor agent (Kiro) is working on this task:\n"
                f"  {original_task}\n\n"
                f"Its recent output:\n{context}\n\n"
                f"Based on this progress, suggest up to 3 concrete next actions. "
                f"If the task appears complete, say DONE.\n"
                f"Format: one action per line, prefixed with ACTION:"
            )
            print(flush=True)
            self._log_event("proposal_requested", task_id, "copilot", cycle=cycle)

            # Extract proposals from streamed text
            proposals = self._extract_proposals()
            if not proposals:
                log.info("No proposals or task complete. Stopping.")
                break

            self._display_proposals(proposals, cycle)
            approved = self._get_approval(proposals)
            if not approved:
                log.info("No proposals approved. Stopping.")
                break

            # Execute approved actions via Kiro
            for action in approved:
                log.info("Dispatching to Kiro: %s", action[:80])
                self.kiro.drain_streamed_text()

                kiro_resp = self.kiro.prompt(action)
                kiro_text = self.kiro.drain_streamed_text()
                self._kiro_output.append(kiro_text)
                print(flush=True)

                self._log_event("proposal_dispatched", task_id, "kiro",
                                content=action[:100])

                if "error" in kiro_resp:
                    log.warning("Kiro action failed: %s — continuing with remaining actions",
                                kiro_resp.get("error"))

    def _extract_proposals(self) -> list[str]:
        """Parse ACTION: lines from Copilot's streamed ACP output."""
        text = self.copilot.drain_streamed_text()
        if not text:
            return []

        # Check if Copilot says task is done — scan the full response, not just
        # the last 3 lines, so DONE is never missed when buried mid-response.
        if "DONE" in text.upper():
            return []

        actions = []
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.upper().startswith("ACTION:"):
                actions.append(stripped[7:].strip())

        # Fallback: lines starting with $ or >
        if not actions:
            for line in text.splitlines():
                stripped = line.strip()
                if stripped.startswith("$ "):
                    actions.append(stripped[2:])
                elif stripped.startswith("> ") and not stripped.startswith("> ["):
                    actions.append(stripped[2:])
        return actions

    def _display_proposals(self, proposals: list[str], cycle: int):
        print(f"\n{'─'*60}")
        print(f"  Copilot proposals (cycle {cycle}):")
        print(f"{'─'*60}")
        for i, action in enumerate(proposals, 1):
            print(f"  [{i}] {action}")
        print()

    def _get_approval(self, proposals: list[str]) -> list[str]:
        if self.auto_approve:
            print("  [auto-approve] Sending all to Kiro\n")
            return proposals
        return [proposals[i] for i in self._prompt_indices(len(proposals))]

    def _prompt_indices(self, count: int) -> list[int]:
        try:
            choice = input("  Send which? (a=all, n=none, 1,2,...): ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            return []
        if choice == "a":
            return list(range(count))
        if choice == "n" or not choice:
            return []
        indices = []
        for part in choice.split(","):
            try:
                idx = int(part.strip()) - 1
                if 0 <= idx < count:
                    indices.append(idx)
            except ValueError:
                pass
        return indices

    def _log_event(self, event: str, task_id: str, agent: str, **extra):
        entry = {
            "event": event,
            "task_id": task_id,
            "agent": agent,
            "ts": datetime.now(timezone.utc).isoformat(),
            **extra,
        }
        self.event_log.append(entry)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Multi-agent orchestrator: Copilot CLI (ACP) ↔ Kiro CLI (ACP)"
    )
    parser.add_argument("--task", required=True, help="Task for Kiro to execute")
    parser.add_argument("--workdir", default=".", help="Working directory")
    parser.add_argument("--kiro-agent", default=None,
                        help="Kiro custom agent name (e.g. rust-developer)")
    parser.add_argument("--auto-approve", action="store_true",
                        help="Send proposals to Kiro without human approval")
    parser.add_argument("--max-cycles", type=int, default=10,
                        help="Maximum proposal cycles per task (default: 10)")
    parser.add_argument("--log-file", default=None,
                        help="Event log output file (default: ~/.acolite/logs/harness_events.jsonl)")
    args = parser.parse_args()

    workdir = str(Path(args.workdir).resolve())

    # Default log file to ~/.acolite/logs/ to avoid polluting the repo tree
    if args.log_file:
        log_file = Path(args.log_file)
    else:
        log_dir = Path.home() / ".acolite" / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / "harness_events.jsonl"

    # -- Launch both ACP agents -------------------------------------------
    log.info("Starting ACP agents...")
    kiro_proc = launch_kiro_acp(workdir, agent=args.kiro_agent)
    copilot_proc = launch_copilot_acp(workdir)

    try:
        kiro = init_acp_agent(kiro_proc, "kiro", workdir, agent_type="kiro")
        copilot = init_acp_agent(copilot_proc, "copilot", workdir, agent_type="copilot")
    except RuntimeError as e:
        log.error("%s", e)
        kiro_proc.terminate()
        copilot_proc.terminate()
        sys.exit(1)

    # -- Print banner ------------------------------------------------------
    print(f"""
╔══════════════════════════════════════════════════════════╗
║  Orchestrator: Kiro CLI ↔ Copilot CLI (both ACP)        ║
╠══════════════════════════════════════════════════════════╣
║  Task:     {args.task[:47]:<47s} ║
║  Workdir:  {workdir[:47]:<47s} ║
║  Kiro:     {'ACP stdio (kiro-cli acp)':<47s} ║
║  Copilot:  {'ACP stdio (copilot --acp --stdio)':<47s} ║
║  Agent:    {(args.kiro_agent or 'default'):<47s} ║
║  Approve:  {'auto' if args.auto_approve else 'manual':<47s} ║
║  Cycles:   {str(args.max_cycles):<47s} ║
╚══════════════════════════════════════════════════════════╝
""")

    # -- Run orchestration -------------------------------------------------
    orch = Orchestrator(kiro, copilot, workdir, args.auto_approve,
                        max_cycles=args.max_cycles)

    try:
        final = orch.run(args.task)
    except KeyboardInterrupt:
        log.info("Interrupted — shutting down")
        final = {"status": "canceled", "events": len(orch.event_log)}
    finally:
        # Write event log
        with open(log_file, "w") as f:
            for event in orch.event_log:
                f.write(json.dumps(event) + "\n")
        log.info("Event log: %s (%d events)", log_file, len(orch.event_log))

        # Shutdown both agents
        kiro.shutdown()
        copilot.shutdown()

    # -- Summary -----------------------------------------------------------
    status = final.get("status", "unknown")
    print(f"\n{'─'*58}")
    print(f"  Task status: {status}")
    print(f"  Events:      {final.get('events', len(orch.event_log))}")
    print(f"  Log:         {log_file}")
    print(f"{'─'*58}\n")

    sys.exit(0 if status == "completed" else 1)


if __name__ == "__main__":
    main()
