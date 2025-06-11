%%
figure 
hold on 
semilogy(errStruct_stephen.f)
semilogy(errStruct_nic.f, 'r--')
legend({'stephen', 'nic'})