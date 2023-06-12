function [alpha,alpha_threshold] = alpha_update(alpha,alpha_threshold,theta,blocks_indicator, powell_factor)

if all(alpha == alpha_threshold) && ~any(blocks_indicator)
    alpha = theta*alpha;
    alpha_threshold = powell_factor*alpha_threshold;
end

end

