%%
figure 
hold on 
semilogy(errStruct_stephen.f, 'LineWidth',5)
semilogy(errStruct_nic.f, 'r--','LineWidth',5)
semilogy(errStruct_pq.f, 'k.','LineWidth',5)
legend({'stephen', 'nic', 'proxQN'})