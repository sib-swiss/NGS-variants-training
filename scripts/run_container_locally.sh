docker run \
--rm \
-e JUPYTER_ENABLE_LAB=yes \
-v $PWD:/home/jovyan \
-p 8888:8888 \
-e GRANT_SUDO=yes \
--user root \
geertvangeest/ngs-variants-jupyter:latest \
start-notebook.sh
