%% plot ctilde

if figps_on 
    axes('FontSize',12)                 % change font size of figures
end
policy_stoch=[policy_l 1./policy_u];
plot(policy_stoch(1,:),ppval(c_tilde_spline,(policy_stoch(1,:))),'r','LineWidth',1+figps_on)
hold on
if (rho == ((1-gamma)^(1-alpha*theta)-1)*lambda-(1-alpha*theta)*delta)
    ctilde_true = (1-gamma)^alpha;
    plot(policy_stoch(1,:),repmat(ctilde_true,length(policy_stoch(1,:)),1),'k:')
    axis([0 max(kss,kss_determ)*1.5 -inf .9488])
elseif (alpha == theta)
    ctilde_true = (1-gamma);
    plot(policy_stoch(1,:),repmat(ctilde_true,length(policy_stoch(1,:)),1),'k:')
    axis([0 max(kss,kss_determ)*1.5 -inf .9])
else
    axis([0 max(kss,kss_determ)*1.5 .915 inf])
end
plot(kss,ppval(c_tilde_spline,kss),'rx','MarkerSize',8)
ylabel('ctilde')
xlabel('capital')
hold off