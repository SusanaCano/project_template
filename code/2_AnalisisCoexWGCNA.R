##### ANALISIS DE COEXPRESIÓN USANDO EL PAQUETE WGCNA #####

#Usamos la funcion WGCNA para elegir los genes de mayor correlación

WGCNA.res2 <- WGCNA(exprs.1=expresion3, exprs.2=expresion4, power = 12, variant = "WGCNA")
WGCNA.res2[1:24]
Links2 <- qLinkfilter(expresion3, expresion4, 0.25)
names(Links2)
umbral_ex3<-Links2$rth.1
#vemos que tiene una correlacion de 0.994
umbral_ex4<-Links2$rth.2
#vemos que tiene una correlacion de 0.997


WGCNA.res3 <- WGCNA(exprs.1=expresion3, exprs.2=expresion5, power = 12, variant = "WGCNA")
WGCNA.res3[1:24]
Links3 <- qLinkfilter(expresion3, expresion5, 0.25)
names(Links3)
umbral_ex3<-Links3$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex5<-Links3$rth.2
#vemos que tiene una correlacion de 0.9833

WGCNA.res4 <- WGCNA(exprs.1=expresion3, exprs.2=expresion6, power = 12, variant = "WGCNA")
WGCNA.res4[1:24]
Links4 <- qLinkfilter(expresion3, expresion6, 0.25)
names(Links4)
umbral_ex3<-Links4$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex6<-Links4$rth.2
#Vamos que tiene una correlacion de 0.9291


WGCNA.res5 <- WGCNA(exprs.1=expresion3, exprs.2=expresion7, power = 12, variant = "WGCNA")
WGCNA.res5[1:24]
Links5 <- qLinkfilter(expresion3, expresion7, 0.25)
names(Links5)
umbral_ex3<-Links5$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex7<-Links5$rth.2
#Vamos que tiene una correlacion de 0.990


WGCNA.res6 <- WGCNA(exprs.1=expresion3, exprs.2=expresion8, power = 12, variant = "WGCNA")
WGCNA.res6[1:24]
Links6 <- qLinkfilter(expresion3, expresion8, 0.25)
names(Links6)
umbral_ex3<-Links6$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex8<-Links6$rth.2
#Vamos que tiene una correlacion de 0.996

WGCNA.res7 <- WGCNA(exprs.1=expresion3, exprs.2=expresion9, power = 12, variant = "WGCNA")
WGCNA.res7[1:24]
Links7 <- qLinkfilter(expresion3, expresion9, 0.25)
names(Links7)
umbral_ex3<-Links7$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex9<-Links7$rth.2
#Vamos que tiene una correlacion de 0.0569

#no ponemos la expresion 9 porque tiene una correlacion de 0.05
datos_expresiones<-rbind(expresion3,expresion4,expresion7,expresion8)


