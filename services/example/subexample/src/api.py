"""
@Goal: Provide an example API for STARK CLI modules
@Author: Samuel Nicaise
@Date: May 2026

For the development of a new STARK module, this is what you have to change in this file:
- replace the subexample_main import by your actual module main
- replace RunRequest parameters so they match the parameters of your main function. Keep the generic API parameters (module_name, threads, memory).
- replace the content of the coroutine variable in the run_endpoint function by a call to your main function, with the parameters from the request


This example API can be called 
# From inside a container in the same docker network
curl -X POST --noproxy stark-module-example-submodule-subexample-service-cli "http://stark-module-example-submodule-subexample-service-cli:8000/run" -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/tmp/run","genome":"hg38"}'

# From the local server
curl --noproxy localhost http://localhost:9999/run -H "Content-Type: application/json" -H "Authorization: Bearer key" -d '{"run_dir":"/tmp/run","genome":"hg38"}'

Where 
8000 is the internal port 
9999 the external port
key the API key
See the STARK.docker-compose.yml and STARK.env for the actual values of these parameters.
"""

import asyncio
import os

from fastapi import FastAPI, HTTPException, Header
from pydantic import BaseModel
import uvicorn

from subexample import main as subexample_main

app = FastAPI()
tasks = {}  # used to store running tasks, to be able to cancel them


class RunRequest(BaseModel):
    run_dir: str
    genome: str
    verbosity: str = "info"


class CancelRequest(BaseModel):
    task_id: str


def verify_api_key(authorization_header: str) -> None:
    if not authorization_header.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Invalid Authorization header")

    api_key = authorization_header.split("Bearer ")[1]

    expected_key = os.getenv("SUBMODULE_API_KEY")
    if expected_key is None:
        raise HTTPException(status_code=500, detail="API key not configured")

    if api_key != expected_key:
        raise HTTPException(status_code=403, detail="Invalid API key")


# See this = Header(...) in the function arguments? It means that this parameter is expected to be in the header of the POST request, and it is required (if it was Header(None) it would be optional). FastAPI automatically parses the header and pass it to the function.
@app.post("/run")
async def run_endpoint(request: RunRequest, authorization: str = Header(...)) -> dict:
    verify_api_key(authorization)

    task_id = f"task-{len(tasks) + 1}"

    func = subexample_main
    args = (request.run_dir, request.genome, request.verbosity)
    coroutine = run_task_with_cancellation(func, args)

    # With asyncio, creating a task automatically puts it in the event loop, so it starts running immediately.
    task = asyncio.create_task(coroutine)
    tasks[task_id] = task
    return {"task_id": task_id, "status": "started"}

async def run_task_with_cancellation(func: callable, args: tuple):
    """
    This runs the function subexample_main in a separate thread, so that it doesn't block the API.
    subexample_main arguments are directly passed in the call to to_thread
    """
    try:
        await asyncio.to_thread(func, *args)
    except asyncio.CancelledError:
        # Handle cleanup if necessary
        print("Task was cancelled")
        raise

@app.post("/cancel")
async def cancel_task(request: CancelRequest, authorization: str = Header(...)) -> dict:
    verify_api_key(authorization)

    task_id = request.task_id
    task = tasks.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    if task.done():
        return {"task_id": task_id, "status": "already completed"}
    # When we cancel a task in asyncio, it eventually raises a CancelledError in the task, which stops it. It isn't immediate so we await the task to ensure it has properly stopped.
    task.cancel()
    try:
        await task
    except asyncio.CancelledError:  # NOSONAR
        # It is normal to not raise anything after this error is caught.
        # If there is cleanup to do after cancelling, it can be done here.
        pass
    except Exception as e:
        raise RuntimeError(f"Task cancellation failed: {e}")
    del tasks[task_id]
    return {"task_id": task_id, "status": "cancelled"}


@app.get("/tasks")
async def list_tasks(authorization: str = Header(...)) -> dict:
    """
    The intended use of this endpoint is to be able to manually get task ids to be able to cancel them with the /cancel endpoint, in case you didn't store the task id when you called the /run endpoint.
    """
    verify_api_key(authorization)

    if tasks == {}:
        return {"running_tasks": {}, "completed_tasks": {}}

    running_tasks = {
        task_id: "running" for task_id, task in tasks.items() if not task.done()
    }
    completed_tasks = {
        task_id: "completed" for task_id, task in tasks.items() if task.done()
    }
    return {"running_tasks": running_tasks, "completed_tasks": completed_tasks}


if __name__ == "__main__":
    # While developping, you can add reload=True to automatically restart the server when you edit the code
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)
