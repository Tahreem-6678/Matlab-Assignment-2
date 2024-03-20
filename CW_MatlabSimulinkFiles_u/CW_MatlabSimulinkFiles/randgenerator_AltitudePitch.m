%Only run this program once and save its results to void confusion
clear all
%Generated altitude is a random number between 620 km and 680 km
Altrange=620 + (680-620)*rand(10,1);
Altituderand=round(Altrange(1))
Pitchrange=15 + (30-15)*rand(10,1);
Pitchrand=round(Pitchrange(1))

save myMat.mat Altituderand Pitchrand