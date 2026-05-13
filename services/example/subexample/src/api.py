"""
Example API
This can be called with 
curl -X POST "http://stark-module-example-submodule-subexample-service-cli:8000/run" -H "Content-Type: application/json" -d '{"run_dir":"/tmp/run","genome":"hg38"}'

From Python this could be called with
import httpx
resp = httpx.post("http://stark-module-example-submodule-subexample-service-cli:8000/run", json={"run_dir":"/tmp/run","genome":"hg38"})
print(resp.status_code, resp.json())

If relevant, replace port 8000 by the exterior port that is defined for the container
"""

from fastapi import FastAPI
from pydantic import BaseModel
import uvicorn

from subexample import main as subexample_main

app = FastAPI()

class RunRequest(BaseModel):
    run_dir: str
    genome: str


@app.post("/run")
async def run_endpoint(request: RunRequest) -> dict:
    try:
        subexample_main(request.run_dir, request.genome)
        result = "success"
    except Exception as e:
        result = f"internal module error: {e}"
    return {"result": result}

if __name__ == "__main__":
    # While developping, you can add reload=True to automatically restart the server when you edit the code
    uvicorn.run("main:app", host="0.0.0.0", port=8000)