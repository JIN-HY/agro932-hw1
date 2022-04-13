cc <- read.csv("data/ppj220030-sup-0002-tables1.csv")
table(cc$date)
### add replication information
cc$Rep <- "Rep2"
cc[cc$Row< 3000,] $Rep <- "Rep1"
s1 <- subset(cc, date %in% "Sept1")
fit <- lm(Canopy_Coverage ~ Genotype + Treatment + Genotype:Treatment + Rep, 
          data=s1)
summary(aov(fit))

date.h2 <- function(Date){
  subcc <- subset(cc, date %in% Date)
  fit <- lm(Canopy_Coverage ~ Genotype + Treatment + Genotype:Treatment + Rep, 
            data=subcc)
  print(summary(aov(fit)))
  MS <- summary(aov(fit))[[1]]$'Mean Sq'
  MSp <- MS[1]
  MSpe <- MS[4]
  VA <- (MSp-MSpe)/4
  h2 <- VA/(VA+MSpe/4)
}

dates <- c("July6","Aug12","Aug14","Aug16","Aug20","Aug22","Aug23","Aug26","Aug30","Sept1","Sept3","Sept5")
date_h2 <- sapply(dates, date.h2)
pdf("graphs/hw2h2.pdf", width = 15)
barplot(date_h2, ylab = "heritability", xlab = "date", main = "Narrow sense heritability of canopy coverage over time", ylim = c(0,1))
dev.off()
