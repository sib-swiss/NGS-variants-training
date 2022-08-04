docker run \
--rm \
-p 8443:8443 \
-e PUID=1000 \
-e PGID=1000 \
-e DEFAULT_WORKSPACE=/config/workdir \
-v $PWD:/config/workdir \
geertvangeest/ngs-variants-vscode
