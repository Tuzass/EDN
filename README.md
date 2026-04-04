# EDN - Equações Diferenciais Numéricas

## TP1

Análise de métodos numéricos para resolução de Equações Diferenciais Ordinárias de primeira ordem. Foram implementadas variações do método de Runge-Kutta de 2 estágios (Heun e ponto médio), e o método de Runge-Kutta de 4 estágios (RK4). Para visualizar os resultados dos métodos, deve-se compilar o arquivo `calculator.cpp` e então executar o script `main.py`. No Windows,

``` bash
g++ -O3 src\calculator.cpp -o calculator
python src\main.py
```

No Linux, esses comandos são:
``` bash
g++ -O3 src/calculator.cpp -o calculator
python3 src/main.py
```

O arquivo `settings.ini` contém opções de configuração que permitem a escolha do método, intervalo, tamanho do passo, etc. Note que a leitura desse arquivo é extremamente simples, então valores inadequados muito provavelmente resultarão em erros. Durante a execução do script `main.py`, o programa `calculator.exe` é executado e cria 3 arquivos temporários contendo os dados relacionados à execução do método escolhido nos parâmetros definidos. Após a exibição dos gráficos por meio do pacote _matplotlib_, esses arquivos são removidos automaticamente pelo script. Note também que só é possível utilizar um método por vez, e que não há necessidade de recompilação do arquivo `calculator.cpp` a menos que a função a ser analisada mude ou um novo método seja adicionado.
