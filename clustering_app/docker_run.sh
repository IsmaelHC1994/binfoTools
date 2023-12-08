# docker build -t ihenarejos/clustering-from-genotypes --no-cache .
docker run -it -d --name dendrograms -v /home/ihenarejos/workspace/phd_summary/docker_make/clustering_from_genotypes/results/ ihenarejos/clustering-from-genotypes
docker exec -it dendrograms /bin/bash
