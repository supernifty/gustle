go build *.go
cp gustle gustle_osx
GOOS=linux go build *.go
cp gustle gustle_linux
