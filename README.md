# HDR_flask

## run server
python3 -m flask run

## make docker container
docker build -t container/dockerfile:vflask1 ./container/
docker run -p 5000:5000 -d container/dockerfile:vflask1
docker cp /path/to/your/local/file.txt container_name:/path/in/container/
docker run -p 5000:5000 -v ${PWD}\outputs_flask:/src/static/outputs -d dockerfile:vflask1 
docker build -t dockerfile:vflask1 ./    


## Check code quality
python3 -m black main.py
python3 lint.py --threshold 7

## start unit-tests
python3 -m pytest unit_tests.py