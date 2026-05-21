# Hypergeometric Sums Lab (WZ / Zeilberger-style)

Objetivo: transformar conjecturas numericas em identidades com certificado por recorrencia.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Familias disponiveis

- `binom-sum`: \(\sum_{k=0}^{n} \binom{n}{k} = 2^n\)
- `binom-square`: \(\sum_{k=0}^{n} \binom{n}{k}^2 = \binom{2n}{n}\)
- `weighted-binom`: \(\sum_{k=0}^{n} k\binom{n}{k} = n2^{n-1}\)
- `franel-cube`: \(\sum_{k=0}^{n} \binom{n}{k}^3\) (Franel, recorrencia de ordem 2)

## 1) Certificado WZ para uma familia

Executar:

```bash
python3 src/wz_lab.py prove-family --family binom-square
```

O script:

- busca automaticamente um certificado racional `R(n,k)` em um espaco de ansatz,
- verifica simbolicamente o residual da equacao WZ.

Para `franel-cube`, o buscador de certificado atual (modo normalizado de primeira ordem) retorna `n/a`.
Agora o `prove-family` tenta automaticamente um modo de telescoping de ordem maior usando a recorrencia conhecida da familia.

## 2) Conjectura -> recorrencia

Executar:

```bash
python3 src/wz_lab.py guess-family --family franel-cube --n-start 0 --n-end 30 --order -1 --degree 2
```

Esse comando:

- gera termos exatos da soma,
- tenta adivinhar recorrencia linear com coeficientes polinomiais,
- valida em termos exatos adicionais.

Obs: `--order -1` usa a ordem recomendada da familia.

## 3) Relatorio automatico (3 familias)

Executar:

```bash
python3 src/wz_lab.py report --n-end 30 --degree 2 --max-degree 4
```

Saida:

- status da recorrencia conjecturada,
- ordem de recorrencia usada por familia,
- status do certificado WZ encontrado,
- verificacao simbolica do residual.

## Observacao

- O comando `prove-family` gera um certificado simbolico no estilo WZ para a familia escolhida.
- O comando `guess-family` e util para descoberta de recorrencias candidatas (pipeline de conjectura).
- O comando `report` automatiza o ciclo em todas as familias internas.
- Use `--no-auto-order` no `report` para forcar uma ordem unica em todas as familias.
- Para familias com recorrencia de ordem > 1 (ex.: `franel-cube`), o sistema tenta certificado de telescoping de ordem maior via ansatz racional.
