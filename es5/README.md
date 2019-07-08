# Monte_Carlo exponential integral Esercizio 5.

Horizontal Rule
***

Il seguente programma effettua l'integrazione della funzione esponenziale 
tra 0 e 1. Vengono implementati anche diversi metodi per la riduzione
della varianza degli estimatori monte carlo: 
1. Stratified Sampling
2. Important Sampling
3. Antithetic variates method 

Horizontal Rule
***

Il programma si avvale dell'uso di tre librerie:

1. Random, genera numeri casuali usando la classe TRandom3 di root che
ha un periodo sufficientemente alto per essere considerato un buon 
generatore di numeri casuali. I metodi implementati permettono di generare
numeri numeri casuali in [0,1], numeri stratificati in [0,1] secondo un intervallo 
oppure numeri numeri distribuiti secondo un esponenziale e 2.5*x^1.5 per 
l'important sampling.

2. Function, supporta la valutazione delle funzioni di interesse e 
computa l'integrale secondo monte carlo.

3. Monte_carlo, parla con le librerie Random e Function per avviare la simulazione.
I metodi permettono di simulare vari monte carlo con tecniche di riduzione della varianza
citati in precedenza.


Il Main è il programma principale per l'analisi che genera i plot di interesse per visualizzare
i risultati della simulazione.

Per avviare il programma:
'./compile.sh'
'bash compile.sh'

L'executable sarà situato in '  ./Monte_Carlo/build/main' e sarà eseguito automaticamente dallo script bash.
