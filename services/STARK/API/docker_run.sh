
# to leuch in main STARK folder


# Source env
source service/.env


# Build docker images Centos
docker build --no-cache -t stark-server-centos . -f service/server/centos.Dockerfile

# Build docker images Alpine
docker build --no-cache -t stark-server . -f service/server/Dockerfile


# docker Centos (for dev)
docker run -p $DOCKER_STARK_PORT:8000 -v /var/run/docker.sock:/var/run/docker.sock -v $DOCKER_STARK_FOLDER_ANALYSES:/STARK/analyses --env-file service/.env -v $(pwd)/service/server/simpleapp.py:/src/simpleapp.py --name=stark-server-centos -ti stark-server-centos bash
# docker Centos (for prod)
docker run -p $DOCKER_STARK_PORT:8000 -v /var/run/docker.sock:/var/run/docker.sock -v $DOCKER_STARK_FOLDER_ANALYSES:/STARK/analyses --env-file service/.env --name=stark-server-centos stark-server-centos

# docker Alpine (for dev)
docker run -p $DOCKER_STARK_PORT:8000 -v /var/run/docker.sock:/var/run/docker.sock -v $DOCKER_STARK_FOLDER_ANALYSES:/STARK/analyses --env-file service/.env -v $(pwd)/service/server/simpleapp.py:/src/simpleapp.py --name=stark-server -ti stark-server ash
# docker Alpine (for prod)
docker run -p $DOCKER_STARK_PORT:8000 -v /var/run/docker.sock:/var/run/docker.sock -v $DOCKER_STARK_FOLDER_ANALYSES:/STARK/analyses --env-file service/.env --name=stark-server stark-server


# Exec Server Centos
docker exec -ti stark-server-centos bash

# Exec Server Alpine
docker exec -ti stark-server bash


# launch server Flask within the docker server
python simpleapp.py -p 8000


# Request for a STARK analysis
curl --request POST --data '{"run":"RUN_TEST_TAG:RUN_TEST_TAG9","sample_filter":"P1408"}' http://0.0.0.0:$DOCKER_STARK_PORT/analysis
curl --request POST --data '{"run":"RUN_TEST_TAG:RUN_TEST_TAG9","sample_filter":"P1338"}' http://0.0.0.0:$DOCKER_STARK_PORT/analysis

# Request for queue
curl --request POST --data '' http://0.0.0.0:$DOCKER_STARK_PORT/queue
http://localhost:$DOCKER_STARK_PORT/queue?format=html

# Request for queue log view
curl --request POST --data '?log=0' http://0.0.0.0:$DOCKER_STARK_PORT/queue

# Request for queue info view
curl --request POST --data '?info=0' http://0.0.0.0:$DOCKER_STARK_PORT/queue
