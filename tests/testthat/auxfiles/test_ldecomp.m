% MATLAB code to generate data for testing LDECOMP class methods
%

clear
clc

% Possible values
% algorithms: "chisq", "ddmoments", "ddrobust"
% dataname: "full", "excluded"

get_results("jm", "full")

function get_results(algorithm, dataname)

   % ranges for significant levels
   alpha = [0.05, 0.10];
   gamma = [0.01, 0.05];

   % prepare data
   X = xlsread("People.xlsx", 1, "B2:M33");

   if dataname == "excluded"
      X([1, 10, 20], :) = [];
      X(:, [3, 12]) = [];
   end

   X = bsxfun(@minus, X, mean(X));
   X = bsxfun(@rdivide, X, std(X));

   fprintf("\n\n++++++++++++++++++++++++++++++++++++++++ \n")
   fprintf("\n%s - %s data\n", algorithm, dataname)
   fprintf("------------------\n")
   for sl = 1:2
      fprintf("\n\nAlpha = %.2f gamma = %.2f\n", alpha(sl), gamma(sl))
      fprintf("*****************************\n")
      [Qlim, T2lim, extremes, outliers] = get_lim(X, algorithm, alpha(sl), gamma(sl));

      fprintf("\nQ limits: alpha = %.2f gamma = %.2f\n", alpha(sl), gamma(sl))
      fprintf("%.8f, ", Qlim(1, :)); fprintf("\n")
      fprintf("%.8f, ", Qlim(2, :)); fprintf("\n")

      fprintf("\nT2 limits:\n")
      fprintf("%.8f, ", T2lim(1, :)); fprintf("\n")
      fprintf("%.8f, ", T2lim(2, :)); fprintf("\n")

      fprintf("\nOutliers (A = 1, 3, 5):\n")
      fprintf("%d, ", find(outliers(:, 1))); fprintf("\n")
      fprintf("%d, ", find(outliers(:, 3))); fprintf("\n")
      fprintf("%d, ", find(outliers(:, 5))); fprintf("\n")

      fprintf("\nExtremes (A = 1, 3, 5):\n")
      fprintf("%d, ", find(extremes(:, 1))); fprintf("\n")
      fprintf("%d, ", find(extremes(:, 3))); fprintf("\n")
      fprintf("%d, ", find(extremes(:, 5))); fprintf("\n")
   end
end


function [Qlim, T2lim, extremes, outliers] = get_lim_pls(X, alpha, gamma, algorithm)

   ncomp = size(X, 2) - 1;

   [~, ~, P] = svd(X);
   T = X * P;
   U = bsxfun(@rdivide, T, sqrt(sum(T.^2) / (size(T, 1) - 1)));

   Qlim = zeros(2, ncomp);
   T2lim = zeros(2, ncomp);
   outliers = zeros(size(X, 1), ncomp);
   extremes = zeros(size(X, 1), ncomp);

   for A = 1:ncomp
      E = X - T(:, 1:A) * P(:, 1:A)';
      Q = sum(E.^2, 2);
      T2 = sum(U(:, 1:A).^2, 2);

      options = struct();
      options.algorithm = algorithm;

      Qlim(1, A) = residuallimit(E, 1 - alpha, options);
      Qlim(2, A) = residuallimit(E, (1 - gamma)^(1/numel(Q)), options);

      T2lim(1, A) = tsqlim(numel(Q), A, 1 - alpha);
      T2lim(2, A) = tsqlim(numel(Q), A, (1 - gamma)^(1/numel(Q)));

      outliers(:, A) = (Q >= Qlim(2, A) | T2 >= T2lim(2, A));
      extremes(:, A) = (Q >= Qlim(1, A) | T2 >= T2lim(1, A)) & ~outliers(:, A);
   end
end

function [Qlim, T2lim, extremes, outliers] = get_lim_dd(X, alpha, gamma, robust)
   ncomp = size(X, 2) - 1;

   Qlim = zeros(2, ncomp);
   T2lim = zeros(2, ncomp);
   outliers = zeros(size(X, 1), ncomp);
   extremes = zeros(size(X, 1), ncomp);

   for A = 1:ncomp
      m = DDSimca(X, A);
      m.Alpha = alpha;
      m.Gamma = gamma;
      m.Scaling = 1;
      m.Centering = 1;
      m.BorderType = "chi-square";

      if robust
         m.EstimationMethod = "robust";
      end

      q0(A) = m.OD_mean * (size(X, 2) - 1);
      h0(A) = m.SD_mean * (size(X, 1) - 1);
      Nq(A) = m.DoF_OD;
      Nh(A) = m.DoF_SD;

      Qlim(1, A) = q0(A)  * m.CriticalLevel(2);
      Qlim(2, A) = q0(A)  * m.OutlierLevel(2);
      T2lim(1, A) = h0(A) * m.CriticalLevel(1);
      T2lim(2, A) = h0(A) * m.OutlierLevel(1);
      outliers(:, A) = m.OutlierObjects';
      extremes(:, A) = m.ExtremeObjects';
   end
end

function [Qlim, T2lim, extremes, outliers] = get_lim(X, algorithm, alpha, gamma)
   if algorithm == "chi2" || algorithm == "jm"
      [Qlim, T2lim, extremes, outliers] = get_lim_pls(X, alpha, gamma, algorithm);
   else
      [Qlim, T2lim, extremes, outliers] = get_lim_dd(X, alpha, gamma, (algorithm == "ddrobust"));
   end
end