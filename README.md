LSN_exercises_delivery

Nelle cartelle delle singole lezioni possono comparire le seguenti librerie:
    Vettore:	libreria personale per la creazione e l'utilizzo di vettori
    random:	per la generazione di numeri casuali, arricchita a partire da quella fornita sul sito Ariel
    Experiment:	contiene oggetti e metodi di base (alcuni virtuali) su cui si basano le successive librerie, più specifiche
    MolDyn_NVE:	libreria a oggetti sviluppata a partire dal codice fornito per simulazioni di MD
    Metropolis:	derivata da Experiment, contiene oggetti e metodi generici per analizzare sistemi sfruttando l'algoritmo di Metropolis
    Hydrogen:	derivata da Metropolis, libreria a oggetti per la misura del raggio medio di un atomo di H



NOTA sull'implementazione del blocking method:
Fino alla lezione03 per il blocking method progressivo salvo prima tutti i risultati dei singoli blocchi in un vettore di N componenti e poi ne calcolo media e deviazione standard della media salvandoli in altri due vettori (le cui componenti descrivono quindi l'evoluzione di questi valori all'aumentare del numero di blocchi utilizzati).
Questo non è necessario e per risparmiare spazio in memoria si possono invece definire delle variabili cumulative su cui poi calcolare le grandezze di interesse e stampare queste ultime direttamente in un file di output, così come faccio nelle lezioni successive.
Per queste lezioni iniziali ho però ritenuto non necessario implementare questo accorgimento: la quantità di dati in esse analizzata è sempre relativamente ridotta ed è trascurabile la possibilità sia di avere problemi di memoria sia di operare in questo senso per migliorare leprestazioni.
In questa prima fase del laboratorio ho quindi preferito implementare un metodo che fosse immediato e comodo da chiamare a livello di main come semplice metodo della classe DataVett piuttosto che concentrarmi sul risparmio di memoria.



README lezione04:
Per produrre i risultati necessari al jupyter notebook basta lanciare main02.sh e main03.sh, che raccolgono una serie di comandi bash per la corretta esecuzione in sequenza dell'insieme di ripetizioni dei tre main necessari a soddisfare le richieste degli esercizi. 
I file nella cartella frame sono invece stati prodotti con lanci singoli (start.sh e restart.sh sono utili per avere le corrette condizioni iniziali o di ripartenza) del main01.x (nstep=1000) e sono utili per rappresentazioni dinamiche con ovito.
