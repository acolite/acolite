//! Secure credential management
//!
//! Sources (in priority order):
//! 1. Environment variables: EARTHDATA_u, EARTHDATA_p
//! 2. .netrc file (Python ACOLITE standard)
//! 3. Config file (~/.acolite/credentials.txt)

use crate::{Result, AcoliteError};
use std::path::PathBuf;
use zeroize::Zeroize;

/// Credential source
#[derive(Debug, Clone, PartialEq)]
pub enum CredentialSource {
    EnvVars,
    Netrc,
    ConfigFile(PathBuf),
}

/// Secure credential container — password zeroized on drop
pub struct Credentials {
    pub username: String,
    password: String,
    token: Option<String>,
    pub source: CredentialSource,
}

impl Drop for Credentials {
    fn drop(&mut self) {
        self.password.zeroize();
        if let Some(ref mut t) = self.token {
            t.zeroize();
        }
    }
}

impl Credentials {
    pub fn load() -> Result<Self> {
        Self::from_env()
            .or_else(|_| Self::from_netrc())
            .or_else(|_| Self::from_config())
            .map_err(|_| AcoliteError::Config(
                "No EarthData credentials found. Set EARTHDATA_u/EARTHDATA_p, \
                 add 'machine earthdata' to ~/.netrc, or create ~/.acolite/credentials.txt".into()
            ))
    }

    pub fn password(&self) -> &str { &self.password }
    pub fn token(&self) -> Option<&str> { self.token.as_deref() }

    fn from_env() -> Result<Self> {
        Ok(Self {
            username: std::env::var("EARTHDATA_u").map_err(|_| AcoliteError::Config("".into()))?,
            password: std::env::var("EARTHDATA_p").map_err(|_| AcoliteError::Config("".into()))?,
            token: std::env::var("EARTHDATA_TOKEN").ok(),
            source: CredentialSource::EnvVars,
        })
    }

    fn from_netrc() -> Result<Self> {
        let path = dirs::home_dir()
            .ok_or_else(|| AcoliteError::Config("no home".into()))?
            .join(".netrc");
        let content = std::fs::read_to_string(&path)
            .map_err(|_| AcoliteError::Config("no .netrc".into()))?;

        let mut found = false;
        let (mut user, mut pass) = (None, None);
        let mut tokens = content.split_whitespace();

        while let Some(tok) = tokens.next() {
            match tok {
                "machine" => found = tokens.next() == Some("earthdata"),
                "login" if found => user = tokens.next().map(String::from),
                "password" if found => pass = tokens.next().map(String::from),
                _ => {}
            }
        }

        Ok(Self {
            username: user.ok_or_else(|| AcoliteError::Config("no login".into()))?,
            password: pass.ok_or_else(|| AcoliteError::Config("no password".into()))?,
            token: None,
            source: CredentialSource::Netrc,
        })
    }

    fn from_config() -> Result<Self> {
        for path in [
            dirs::home_dir().map(|h| h.join(".acolite/credentials.txt")),
            Some(PathBuf::from("config/credentials.txt")),
            dirs::home_dir().map(|h| h.join(".easi-workflows-auth.conf")),
        ].into_iter().flatten() {
            if let Ok(c) = Self::parse_ini(&path) { return Ok(c); }
        }
        Err(AcoliteError::Config("no config file".into()))
    }

    fn parse_ini(path: &PathBuf) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .map_err(|_| AcoliteError::Config("unreadable".into()))?;
        let (mut user, mut pass, mut tok, mut in_sec) = (None, None, None, false);
        for line in content.lines() {
            let line = line.trim();
            if line == "[earthdata]" { in_sec = true; continue; }
            if line.starts_with('[') { in_sec = false; continue; }
            if in_sec {
                if let Some((k, v)) = line.split_once(':') {
                    match k.trim() {
                        "user" | "login" => user = Some(v.trim().to_string()),
                        "password" => pass = Some(v.trim().to_string()),
                        "token" => tok = Some(v.trim().to_string()),
                        _ => {}
                    }
                }
            }
        }
        Ok(Self {
            username: user.ok_or_else(|| AcoliteError::Config("no user".into()))?,
            password: pass.ok_or_else(|| AcoliteError::Config("no password".into()))?,
            token: tok,
            source: CredentialSource::ConfigFile(path.clone()),
        })
    }
}

impl std::fmt::Debug for Credentials {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Credentials")
            .field("username", &self.username)
            .field("password", &"[REDACTED]")
            .field("token", &self.token.as_ref().map(|_| "[REDACTED]"))
            .field("source", &self.source)
            .finish()
    }
}

/// AWS profile from environment only
pub fn aws_profile() -> Option<String> {
    std::env::var("AWS_PROFILE").ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debug_redacts() {
        let c = Credentials {
            username: "user".into(), password: "s3cretpass".into(),
            token: Some("bearer_xyz".into()), source: CredentialSource::EnvVars,
        };
        let s = format!("{:?}", c);
        assert!(!s.contains("s3cretpass"));
        assert!(!s.contains("bearer_xyz"));
        assert!(s.contains("[REDACTED]"));
    }
}
