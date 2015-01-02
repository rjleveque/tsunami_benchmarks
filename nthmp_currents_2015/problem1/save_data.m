
% convert data to ascii text files for reading into Python:

load ../all_data/1_Submerged_Island_Lab/comparison_data/SL_S1_U.dat;
te1 = SL_S1_U(:,1);
ue1 = SL_S1_U(:,2);
s1u = [te1, ue1];
save('S1u.txt','s1u','-ascii')

load ../all_data/1_Submerged_Island_Lab/comparison_data/SL_S1_V.dat;
te2 = SL_S1_V(:,1);
ve2 = SL_S1_V(:,2);
s1v = [te2, ve2];
save('S1v.txt','s1v','-ascii')

load ../all_data/1_Submerged_Island_Lab/comparison_data/SL_S2_U.dat;
te3 = SL_S2_U(:,1);
ue3 = SL_S2_U(:,2);
s2u = [te3, ue3];
save('S2u.txt','s2u','-ascii')

load ../all_data/1_Submerged_Island_Lab/comparison_data/SL_S2_V.dat;
te4 = SL_S2_V(:,1);
ve4 = SL_S2_V(:,2);
s2v = [te4, ve4];
save('S2v.txt','s2v','-ascii')


