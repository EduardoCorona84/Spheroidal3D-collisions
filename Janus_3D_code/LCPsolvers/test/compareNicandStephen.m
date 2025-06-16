%%

minval = min([errStruct_stephen.f, errStruct_nic.f', errStruct_pq.f']) - 10; 
stephen = errStruct_stephen.f - minval;
nic = errStruct_nic.f - minval;
pq =  errStruct_pq.f - minval;
figure 
hold on 
plot(stephen, 'LineWidth',5)
plot(nic, 'r--','LineWidth',5)
plot(pq, 'k.','LineWidth',5)
set(gca, 'YScale', 'log')
legend({'stephen', 'nic', 'proxQN'})