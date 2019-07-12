%% plot errors of the last iteration
if figps_on 
    axes('FontSize',12)                 % change font size of figures
end
semilogy([policy_l(1,:) 1./policy_u(1,:)],abs([(policy_l(2,:)-policy_compare_last_iteration_l)./policy_compare_last_iteration_l (1./policy_u(2,:)-1./policy_compare_last_iteration_u)./policy_compare_last_iteration_u]),'LineWidth',1+figps_on)
hold on
semilogy([kss-eps kss],[10e-20 10e-4],'r:','LineWidth',1)
if (alpha ~= theta)
    axis([0 max(kss,kss_determ)*1.5 10e-15 10e-4])
else
    axis([0 max(kss,kss_determ)*1.5 10e-18 10e-4])
end
ylabel('relative change last iteration (log scale)')
xlabel('capital')
hold off
