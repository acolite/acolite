## def check task
## checks until/whether Google Drive export task is complete
## written by Quinten Vanhellemont, RBINS
## 2022-04-11

def check_task(task, task_sleep = 10):
    import time
    while task.status()['state'] in ['READY', 'RUNNING']:
        status = task.status()
        print('Task {} state {}: sleeping {}'.format(status['description'], status['state'], task_sleep))
        time.sleep(task_sleep)
    status = task.status()
    print('Task completed in {:.1f} seconds with state {}'.format((status['update_timestamp_ms'] - status['creation_timestamp_ms'])/1000, status['state']))
    return(status)
