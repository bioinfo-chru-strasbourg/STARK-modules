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
curl -X POST --noproxy stark-module-example-submodule-subexample-service-cli http://stark-module-example-submodule-subexample-service-cli:8000/run -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/test/run","fibonacci_n":"100"}'
```

- From the local server

```bash
curl -X POST --noproxy localhost http://localhost:9999/run -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/test/run","fibonacci_n":"100"}'
```

Where:

- 8000 is the internal port
- 9999 is the external port
- "key" is the API key

See the STARK.docker-compose.yml and STARK.env for the actual values of these parameters. Always change the API key in production.

## Other API commands

```bash
# cancel a task
curl -X POST --noproxy localhost http://localhost:9999/cancel -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"task_id": "task-1"}'

# get task list (useful if you forgot to store the id of a task you want to cancel)
curl --noproxy localhost http://localhost:9999/tasks -H "Authorization: Bearer key" 
```
