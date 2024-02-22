docker run \
--rm \
-p 8443:8443 \
-e PUID=1000 \
-e PGID=1000 \
-e DEFAULT_WORKSPACE=/config/project \
-v $PWD:/config/project \
geertvangeest/ngs-variants-vscode
