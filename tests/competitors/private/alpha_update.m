function [alpha,alpha_threshold] = alpha_update(alpha,alpha_threshold,theta,blocks_indicator)

if all(alpha == alpha_threshold) && all(blocks_indicator)
    alpha = theta*alpha;
    alpha_threshold = 0.1*alpha_threshold;
end

end

