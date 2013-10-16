# DANIEL CAMILO OSORIO
# UNIVERSIDAD INDUSTRIAL DE SANTANDER
# d.osorio@me.com
# USO: CN("PROTEINA")

CN<-function(seq){
  require(seqinr)
  seq<-s2c(seq)
  t<-table(seq)
  x<-0
  for(i in 1:length(t)){
    if (names(t[i])=="K"){
      x<-x+as.numeric(t[i])
    }
    if (names(t[i])=="R"){
      x<-x+as.numeric(t[i])
    }
    if (names(t[i])=="D"){
      x<-x-as.numeric(t[i])
    }
    if (names(t[i])=="E"){
      x<-x-as.numeric(t[i])
    }
  }
  x
}