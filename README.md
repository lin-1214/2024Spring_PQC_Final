# Neon Assignment

## Modifications
1. Speed up the polynomial multiplication on `ntruhps2048677`
2. The files we add are placed under `./src/ref/ntruhps2048677` 

## Execution
Follow the below steps to execute the code \
1. Go to the correct directory

```bash
cd ./src/ref/ntruhps2048677
```

2. Compile the file

```bash
make
```

3. Test the correctness of the implementation

```bash
./test
```

4. Test the iteration time of the implementation

```bash
./speed
```

## Reference
* https://eprint.iacr.org/2022/930
* https://github.com/GMUCERG/PQC_NEON
* https://github.com/sixkwnp/PQC_Cortex-M4_NTRU-NTT
* https://github.com/vincentvbh/multi-params-ntt_NTRU_NTRUPrime

