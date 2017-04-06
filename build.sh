go build fastaquery.go
cp fastaquery fastaquery_osx
GOOS=linux go build fastaquery.go
cp fastaquery fastaquery_linux
