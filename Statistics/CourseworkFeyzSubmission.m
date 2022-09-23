%% *Setup*
% *Random seed*
rng(01778525)

% *Import Data*
train_data = readtable("train1778525.csv");
test_data = readtable("test1778525.csv");

% *Store training data in arrays*
y = table2array(train_data(:,1));    % Fuel consumption
mass = table2array(train_data(:,2));
time = table2array(train_data(:,3)); % Acceleration time
disp = table2array(train_data(:,4)); % disp
fuel = table2array(train_data(:,5)); % Engine fuel type
color = table2array(train_data(:,6));

%% *Question 1*
% A) Summary statistics*
y_mean = mean(y)
y_median = median(y)
y_std = std(y)
y_iqr = iqr(y)

% Present in table
headers = ["Mean"; "Median"; "Standard deviation"; "Interquartile Range"];
y_stats = table(y_mean, y_median, y_std, y_iqr)

%% Plot histogram
y_hist = histogram(y,10)
title("Histogram of Fuel consumption (l/100km")

%% *B) i ) Boxplots for categorical variables*

boxplot(y,fuel)
title("Boxplot of Fuel consumption (l/100km) vs Fuel type")
xlabel("Fuel type")
ylabel("Fuel consumption (l/100km)")

boxplot(y,color)
title("Boxplot of Fuel consumption (l/100km) vs Car color")
xlabel("Car color")
ylabel("Fuel consumption (l/100km)")

%% *Encode categorical variable fuel*
fuel = grp2idx(fuel); 

%% *B)ii) Scatterplots of Continuous variables + correlation matrices*

% *Mass*
scatter(mass,y)
xlabel("Car mass(kg)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs mass")
corrcoef(mass,y)
%%
% *Acceleration time*
scatter(time,y)
xlabel("Acceleration time to 100km (s)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs acceleration time")
corrcoef(time,y)

% We cant obviously see a linear relationship, so I will test some transformations
scatter(time.*time,y)
xlabel("Acceleration time to 100km (s^2)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs acceleration time squared")
corrcoef(time.*time,y)

% This does not appear to give a good relationship either
scatter(log(time),y)
xlabel("Acceleration time to 100km (log(s))")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs log acceleration time ")
corrcoef(log(time),y)

%% 
% *Engine displacement*
scatter(disp,y)
xlabel("Engine disp(l)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs Engine disp")
corrcoef(disp,y)

%% *Q2) Modelling Simple linear regression with just mass*
% Start from simple models and work my way up by adding complexity
% *Fit a linear regression for fuel consumption vs mass*

% Fit model and predict
model = fitlm(mass,y)
pred = model.predict(mass)

% Plot
scatter(mass,y)
hold on
plot(mass,pred)
xlabel("Car mass(kg)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs mass")
hold off

%% *B) Evaluate Performance Metrics*
rsquared = model.Rsquared
mse = model.MSE
aic = model.ModelCriterion

%% *C) Test different models*
% First plot some scatter plots of possible interaction terms I can use*
% Mass * fuel
scatter(mass.*fuel,y)
xlabel("Mass * fuel type")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs Mass-FuelType")
corrcoef(mass.*fuel,y)

%% 
% This interaction term shows very strong correlation with fuel consumption
% disp * time
scatter(disp.*time,y)
xlabel("Disp*time")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs Disp*time")
corrcoef(disp.*time,y)

%% 
% Very weak correlation
% disp/time
scatter(disp./time,y)
xlabel("Disp/time")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs Disp/time")
corrcoef(disp./time,y)

%% 
% This shows a slightly higher correlation than disp*time100 so I'll prefer 
% this one over that.
% Mass*disp
scatter(mass.*disp,y)
xlabel("Car mass(kg)*disp(l*kg)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs mass * disp")
corrcoef(mass.*disp,y)

%% 
% Fuel*disp
scatter(fuel.*disp,y)
xlabel("Car mass(kg)*disp(l*kg)")
ylabel("Fuel consumption (l/100km)")
title("Fuel consumption vs mass * disp")
corrcoef(mass.*disp,y)

%% Create table to populate with summary stats for the models I wish to test
sz = [10 5];
varTypes = ["double","string", "double","double","double"];
varNames = ["Model#", "Terms", "R^2", "MSE", "AIC"];
model_summaries = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames)

%% 
% For all below models I first create the x dataset, fit the model, plot it, 
% then add summary stats to the comparison table. 

%% *Model 1) Mass + fuel*
n = 1;
x = [mass,fuel];
model = fitlm(x,y)
hold off
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,fuel",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 2)* Mass, engine disp, and fuel type
n = 2;
x = [mass,disp, fuel];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 3)* Mass, engine disp, fuel, and time100
n = 3;
x = [mass, disp, fuel,time];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,time",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 4)* Mass, engine disp, fuel, and mass*fuel
n = 4;
mass_fuel_inter = mass.*fuel;
x = [mass, disp, fuel ,mass_fuel_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,mass*fuel",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 5)* Mass, engine disp, fuel,time, and mass*fuel
n = 5;
mass_fuel_inter = mass.*fuel;
x = [mass, disp, fuel,time, mass_fuel_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,time,mass*fuel",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 6)* Mass, engine disp, fuel, time, mass*fuel, fuel*disp
n = 6;
mass_fuel_inter = mass.*fuel;
fuel_disp_inter = fuel.*disp;
x = [mass, disp, fuel,time, mass_fuel_inter,fuel_disp_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,time,mass*fuel,fuel*disp",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 7)* Mass, engine disp, fuel, time, mass*fuel, disp/time
n = 7;
mass_fuel_inter = mass.*fuel;
disp_time_inter = disp./time;
x = [mass, disp, fuel,time, mass_fuel_inter,disp_time_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,time,mass*fuel,disp/time",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 8)* Mass, engine disp, fuel, mass*fuel, disp/time
n = 8;
mass_fuel_inter = mass.*fuel;
disp_time_inter = disp./time;
x = [mass, disp, fuel, mass_fuel_inter,disp_time_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,disp,fuel,mass*fuel,disp/time",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% *Model 9)* Engine disp, time, mass*fuel 

n = 9;
mass_fuel_inter = mass.*fuel;
x = [disp, time, mass_fuel_inter];
model = fitlm(x,y)
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "disp,time,mass*fuel",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Model 10)* Mass,fuel, disp,time, + 5 interaction terms generated by matlab
n = 10
x = [mass, fuel, disp,time];
model = fitlm(x,y,"interactions")
plot(model)
title("Model",n)
model_summaries(n,:) = {n, "mass,fuel,diso,time + 5 interactions",model.Rsquared.Ordinary, model.MSE, model.ModelCriterion.AIC};

%% Inspect Model Summaries*
model_summaries

%% *Q3) A)Calculate residuals*
% *Redo model*
mass_fuel_inter = mass.*fuel
x = [mass, disp, fuel,time, mass_fuel_inter];
model = fitlm(x,y)
%%Plot residuals vs fitted values*
plotResiduals(model,"fitted")

%% QQ Plot of Residuals*
residuals = table2array(model.Residuals(:,1))
qqplot(residuals)

%% *Theoretical and Empriical CDF*
plotResiduals(model,"probability")

%% B)Interpret the regression coefficients and suggest improvements*
"X1= mass, X2= disp X3= fuel X4=time X5= mass*fuel"
model

%%  *Intercept*
model_noint = fitlm(x,y,"Intercept",false)
model_noint
model_noint.ModelCriterion.AIC
model_noint.Rsquared

%% *Betas*
mass_fuel_inter = mass.*fuel;
x = [mass, fuel,time, mass_fuel_inter];

model_new = fitlm(x,y)
plot(model_new)
model_new.ModelCriterion.AIC
model_new.Rsquared
model_new.MSE


model = model_new;

%% *Q4) Clustering*
% *Perform cluster analysis to identify two distinct groups of cars*

% Perform K means
[idx,C] = kmeans(x,2);
% Use Logical indexing to separate data by category
x_idx = [idx,x];
cat_1 = x_idx(:,1) < 2 ;
cat_2 = x_idx(:,1) > 1;
x_1 = x(cat_1,:);
x_2 = x(cat_2,:);

%%
% Further indexing two separate each category by fuel type (will be useful)
% in below plot) - 1 is diesel, 2 is petrol
% Create logical indices
idx_1d =x_1(:,2) < 2;
idx_1p = x_1(:,2) > 1;
idx_2d = x_2(:,2) < 2;
idx_2p =x_2(:,2) > 1;
% Separate datasets
x_1_diesel = x_1(idx_1d,:);
x_1_petrol = x_1(idx_1p,:);
x_2_diesel = x_2(idx_2d,:); % No diesel cars in category 2
x_2_petrol = x_2(idx_2p,:);

%% 
% *Plot the different categories on cars in a mass v fuel consumption plot, 
% with fuel types indicated*

% Get the fuel consumptions for each category
y1_diesel = model.predict(x_1_diesel);
y1_petrol = model.predict(x_1_petrol);
y2_diesel = model.predict(x_2_diesel);
y2_petrol = model.predict(x_2_petrol);
% Plot
scatter(x_1_diesel(:,1), y1_diesel, 12, "red","filled")
hold on
scatter(x_1_petrol(:,1), y1_petrol, 12, "magenta","filled")
scatter(x_2_petrol(:,1), y2_petrol, 12, "green","filled")

legend("Category 1-diesel","Category 1-petrol", "Category 2-petrol")
title("Mass vs Fuel consumption")
axis([750,2500,0,10])
hold off

%% Q5) Predictions on test set*

% *Prepare data*
y_test = table2array(test_data(:,1));

% Prepare x terms
x1 = table2array(test_data(:,2)); % Mass
x2 = table2array(test_data(:,5)); % Fuel
x2 = grp2idx(x2); % Encoded fuel
x3 = table2array(test_data(:,3)); % Time
x4 = x1.*x2; %mass *fuel

% Overall x matrix
x_test = [x1,x2,x3,x4];

%% Summary stats for output
y_mean = mean(y_test)
y_median = median(y_test)
y_std = std(y_test)
y_iqr = iqr(y_test)

%% Fit model and plot*
y_pred = predict(model,x_test);
plot(x1,y_test,".",x1,y_pred,"x")
legend("Fitted value", "Real value")
title("Fuel yiciency vs mass(test data)")

%% For comparison we can plot this identical plot for the training data
y_pred_test = predict(model,x);
plot(mass,y_pred_test,".",mass,y,"x")
legend("Fitted value", "Real value")
title("Fuel yiciency vs mass(train data)")

%% Calculate metrics*
mse_test = mean((y_test - y_pred).^2)
mse_train = model.MSE
