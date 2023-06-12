function [alpha,alpha_threshold] = alpha_update(alpha,alpha_threshold,theta)

if all(alpha == alpha_threshold) && 
    alpha = theta*alpha;
    alpha_threshold = 0.1*alpha_threshold;
end

end

