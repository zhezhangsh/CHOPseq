# reverse-complement a PWM matrix 
CHOPseq.ReversePWM<-function(pwm) {
x<-pwm[4:1,];
x<-x[, ncol(x):1];
rownames(x)<-rownames(pwm);
colnames(x)<-colnames(pwm);
x;
} 