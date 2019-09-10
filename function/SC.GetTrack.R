# Retrieve a UCSC track as a data.frame
SC.GetTrack<-function(genome, track) {
# genome, track    Genome and track name
library(rtracklayer);
session <- browserSession();
genome(session) <- genome;
query <- ucscTableQuery(session, track);

getTable(query);
}

SC.GetTrackNames<-function(genome) {
library(rtracklayer);
session <- browserSession();
genome(session) <- genome;
trackNames(session);
}