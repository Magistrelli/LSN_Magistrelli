LSN_exercises_delivery

Nelle cartelle delle singole lezioni possono comparire le seguenti librerie:
    Vettore:	libreria personale per la creazione e l'utilizzo di vettori
    random:	libreria per la generazione di numeri casuali, arricchita a partire da quella fornita sul sito Ariel



NOTA sull'implementazione del blocking method:
Fino alla lezione03 per il blocking method progressivo salvo prima tutti i risultati dei singoli blocchi in un vettore di N componenti e poi ne calcolo media e deviazione standard della media salvandoli in altri due vettori (le cui componenti descrivono quindi l'evoluzione di questi valori all'aumentare del numero di blocchi utilizzati).
Questo non è necessario e per risparmiare spazio in memoria si possono invece definire delle variabili cumulative su cui poi calcolare le grandezze di interesse e stampare queste ultime direttamente in un file di output, così come faccio nelle lezioni successive.
Per queste lezioni iniziali ho però ritenuto non necessario implementare questo accorgimento: la quantità di dati in esse analizzata è sempre relativamente ridotta ed è trascurabile la possibilità sia di avere problemi di memoria sia di operare in questo senso per migliorare leprestazioni.
In questa prima fase del laboratorio ho quindi preferito implementare un metodo che fosse immediato e comodo da chiamare a livello di main come semplice metodo della classe DataVett piuttosto che concentrarmi sul risparmio di memoria.  
