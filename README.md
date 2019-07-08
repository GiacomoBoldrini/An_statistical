# PSEUDO-RANDOM numbers Esercizio 1.

Il primo esercizio propone di generare dei numeri pseudo-causuali distribuiti secondo una distribuzione di landau. Viene quindi chiesto di analizzare test di casualità per i numeri appena generati.
Lo svolgimento parte dalla creazione dei generatori pseudo-casuali prng implementato:
- Linear congruent method
-xorshiro plus 128

Vengono successivamenti testati i generatori grazie alla libreria dieharder che implementa un numero elevato di test di casualità.
Il migliore algoritmo viene quindi utilizzato per generare la Landau attraverso metodo di accept reject.


# 1-, 0- Bosons distributions Esercizio 3.
L'esercizio 3 propone di generare 50000 eventi monte carlo secondo la distribuzione angolare dei bosoni vettori:
$$F(\theta, \phi) = \frac{3}{4\pi}[0.5(1-\alpha) + (0.5)(3\alpha-1)cos(\theta)^2 - \beta sin(\theta)^2cos(2\phi)- \sqrt2\gamma sin(2\theta)cos(\phi)]  $$
Con
$$\alpha = 0.65, \beta = 0.06, \gamma = -0.18 $$
Viene quindi chiesto di fittare la distribuzione con metodo ML unbinned, successivamente di binnare i dati in risoluzione dell'1\% in $\theta, \phi$ e di effettuare un fit $\chi^{2}$.
Nello svolgimento è stato utilizzato il pacchetto RooFit di ROOT sviluppato per l'esperimento BaBar a SLAC. Supporta metodi e oggetti molto interessanti per fittare precisamente le distribuzioni attraverso metodi unbinned ML e chi2, oltrettutto permette di creare qualsiasi pdf, sommare pdf e moltiplicarle, computando automaticamente la normalizzazione e la forma funzionale della pdf.
Viene quindi chiesto di testare l'ipotesi di distribuzione uniforme per bosoni 0- rispetto alla precedente di bosoni 1- sui dati generati dalla seconda. Le ipotesi sono semplici in quanto i parametri sono tutti determinati a priori, per lemma di Neyman Pearson la statistica migliore è data dal rapporto delle likelihood che, come logico, sarà piccolissimo favorendo l'ipotesi 1- rispetto a quella di distribuzione uniforme 0-.



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
