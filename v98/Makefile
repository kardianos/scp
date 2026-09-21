.PHONY: test sizespec h0

test:
	go test ./size/

sizespec:
	go build -o sizespec ./cmd/sizespec

h0: test
	go run ./cmd/sizespec -h0 broad   ../v93/runs/faceC/c_deep2.fcs
	go run ./cmd/sizespec -h0 line6   ../v93/runs/faceA/ringA_strang_live.fcs
	go run ./cmd/sizespec -h0 ceiling /space/scp/v93/fcs/rad_tax_S05.fcs
