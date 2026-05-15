Example STARK module
============

The module is named "example", the submodule is named "subexample".

## Running this service with launcher.py

Everything works as usual.

listener.py -> launcher.py -> subexample.py

## Running this service with an API call

This example API can be called

- From inside a container in the same docker network

```bash
curl -X POST --noproxy stark-module-example-submodule-subexample-service-cli "<http://stark-module-example-submodule-subexample-service-cli:8000/run>" -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/tmp/run","genome":"hg38"}'
```

- From the local server

```bash
curl --noproxy localhost <http://localhost:9999/run> -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/tmp/run","genome":"hg38"}'
```

Where:

- 8000 is the internal port
- 9999 the external port
- key the API key

See the STARK.docker-compose.yml and STARK.env for the actual values of these parameters.

## Other API commands

```bash
# cancel a task
curl --noproxy localhost http://localhost:9999/cancel -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"task_id": "task-1"}'

# get task list (useful if you forgot to store the id of a task you want to cancel)
curl --noproxy localhost http://localhost:9999/tasks -H "Authorization: Bearer key" 
```
