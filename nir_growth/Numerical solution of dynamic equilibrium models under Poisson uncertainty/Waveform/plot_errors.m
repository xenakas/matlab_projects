%% plot errors
if figps_on 
    axes('FontSize',12)                 % change font size of figures
end
if (alpha == theta || rho == ((1-gamma)^(1-alpha*theta)-1)*lambda-(1-alpha*theta)*delta)
    semilogy([policy_l(1,:) 1./policy_u(1,:)],abs([policy_l(2,:)-policy_true_l 1./policy_u(2,:)-1./policy_true_u]),'LineWidth',1+figps_on)
else
    semilogy([policy_l(1,:) 1./policy_u(1,:)],abs([policy_l(2,:)-policy_compare_last_iteration_l 1./policy_u(2,:)-1./policy_compare_last_iteration_u]),'LineStyle','none')
end
hold on
semilogy([kss-eps kss],[10e-20 10e-4],'r:','LineWidth',1)
if (alpha ~= theta)
    axis([0 max(kss,kss_determ)*1.5 10e-11 10e-4])
else
    axis([0 max(kss,kss_determ)*1.5 10e-18 10e-4])
end
ylabel('absolute error (log scale)')
xlabel('capital')
hold off
