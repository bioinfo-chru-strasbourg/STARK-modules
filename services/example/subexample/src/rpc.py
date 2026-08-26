"""
@Goal: Provide an example RPC API for STARK CLI modules
@Author: Samuel Nicaise
@Date: June 2026

For the development of a new STARK module, this is what you have to change in this file:
- replace the subexample_main import by your actual module main
- replace RunRequest parameters so they match the parameters of your main function
- replace the content of the coroutine variable in the run_endpoint function by a call to your main function, with the parameters from the request
"""

import os

import fastapi_jsonrpc as jsonrpc
from pydantic import BaseModel
from fastapi import Body

from fastapi import FastAPI, HTTPException, Header
from pydantic import BaseModel
import uvicorn

from subexample import main as subexample_main


app = FastAPI()


class RunRequest(BaseModel):
    run_dir: str
    fibonacci_n: int
    threads: int = 1
    memory: str = "1G"
    verbosity: str = "info"

def verify_api_key(authorization_header: str) -> None:
    if not authorization_header.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Invalid Authorization header")

    api_key = authorization_header.split("Bearer ")[1]

    expected_key = os.getenv("SUBMODULE_API_KEY")
    if expected_key is None:
        raise HTTPException(status_code=500, detail="API key not configured")

    if api_key != expected_key:
        raise HTTPException(status_code=403, detail="Invalid API key")


@app.post("/run")
async def run_endpoint(request: RunRequest, authorization: str = Header(...)) -> dict:
    """
    Main endpoint to start your module's task.

    See this = Header(...) in the function arguments? It means that this parameter is expected to be in the header of the POST request, and it is required (if it was Header(None) it would be optional). FastAPI automatically parses the header and pass it to the function.
    """
    verify_api_key(authorization)

    subexample_main (
        request.run_dir,
        request.fibonacci_n,
        request.threads,
        request.memory,
        request.verbosity,
    )
    return {"message": "Task started successfully"}




if __name__ == "__main__":
    # While developping, you can add reload=True to automatically restart the server when you edit the code. Remove it in production.
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)
