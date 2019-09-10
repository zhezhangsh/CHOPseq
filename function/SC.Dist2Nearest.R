# Minimal absolute distance to nearest feature
SC.Dist2Nearest<-function(query, features) {
# query, features   GRanges objects
if (length(features)==0 | length(query)==0) {
  rep(NA, length(query));
}
else {
start<-start(query);
end<-end(query);

loc<-unique(c(start(features), end(features)));

d1<-abs(start-loc[1]);
d2<-abs(end-loc[1]);

for (i in 1:length(loc)) {
d1<-pmin(d1, abs(start-loc[i]));
d2<-pmin(d2, abs(end-loc[i]));
}

pmin(d1, d2);
}
}
