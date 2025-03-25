%% plot deterministic versus stochastic policy function

if figps_on 
    axes('FontSize',12)                 % change font size of figures
end
plot(policy_determ(1,:),policy_determ(2,:),'LineStyle','--','LineWidth',1+figps_on,'Color',[0 0 1])
hold on
policy_stoch=[policy_l 1./policy_u];
plot(kss_determ,css_determ,'x','MarkerSize',8,'LineWidth',1)
plot(policy_stoch(1,:),policy_stoch(2,:),'r','LineWidth',1+figps_on)
plot(kss,css,'rx','MarkerSize',8,'LineWidth',1)
if (alpha ~= theta)
    policy_true=(theta-1)/theta*abs(policy_stoch(1,:)).^alpha;
    plot(policy_stoch(1,:),policy_true,'k:')
else
    policy_true = (rho-(exp((1-theta)*log(1-gamma))-1)*lambda+(1-theta)*delta)/theta*abs(policy_stoch(1,:));
    plot(policy_stoch(1,:),policy_true,'k:')
end
axis([0 max(kss,kss_determ)*1.5 0 max(css,css_determ)*1.5])
ylabel('consumption')
xlabel('capital')
hold off
