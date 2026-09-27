# Container usage

Build:

```bash
docker build -t bioinfox/snp2cluster-py:0.1.0 .
```

Run validation:

```bash
docker run --rm -v "$PWD:/work" -w /work bioinfox/snp2cluster-py:0.1.0 validate --config configs/example_config.yaml
```

Run analysis:

```bash
docker run --rm -v "$PWD:/work" -w /work bioinfox/snp2cluster-py:0.1.0 run --config configs/example_config.yaml
```
