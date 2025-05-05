function [FK_aa] = aa_cut(Fk,aa_filter)


%%
% aa_filter = ( K_ >= 1/3*n);

%%
NN = size(aa_filter);
if length(NN) == 2
    FK_aa = aa_filter.*Fk;
else
    FK_aa = aa_filter.*Fk;
end


