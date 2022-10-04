# Builds images
docker build -f base.Dockerfile -t hecate-core .
docker run --name hecate -it hecate-core