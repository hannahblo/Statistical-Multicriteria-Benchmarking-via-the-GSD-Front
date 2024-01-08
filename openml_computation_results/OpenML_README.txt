### Ergebnisse der Berechnungen über Weihnachten 2023


### Was wurde berechnet
Alle Algorithmen unten wurden mit SVM verglichen. 
Als epsilon wurde verwendet:
0. 0
1. 0.25
2. 0.5
3. 0.75
4. 1

Die Daten kommen aus OpenML. Der Filter ist wie folgt:
rel_datasets <- data_all %>% filter(status == "active" &
                                      number.of.classes == 2 &
                                      number.of.features < 100 &
                                      number.of.instances < 1000 &
                                      number.of.instances > 100 &
                                      number.of.instances.with.missing.values == 0 &
                                      max.nominal.att.distinct.values < 5
)



### Algorithen
glmnet = 
kknn = KNN
multinom =
ranger = Random Forest
rpart = CART
xgboost = 


### Die Files im Git
.._computation_time.rds = Berechnungszeit
.._dat_set.rds = exact der Datensatz, der verwendet wurde für die Berechnung
.._result.rds = Ergebnisse der Berechnung (also der einzelnen Permutationen und nicht permutiert)
.._values_teststatistics.jpeg = Verteilung der Teststatistik bassierend auf unterschiedlichen gammas
.._regularization.jpeg = sollte eigentlich der gamma-plot sein --> funktioniert aber leider nicht...


### Rechenzeiten
> readRDS("classif.glmnet_computation_time.rds")
Time difference of 14.43411 hours
> readRDS("classif.kknn_computation_time.rds")
Time difference of 13.51942 hours
> readRDS("classif.multinom_computation_time.rds")
Time difference of 1.38964 days
> readRDS("classif.ranger_computation_time.rds")
Time difference of 3.954684 days
> readRDS("classif.rpart_computation_time.rds")
Time difference of 2.647984 days
> readRDS("classif.xgboost_computation_time.rds")
Time difference of 3.987366 days


