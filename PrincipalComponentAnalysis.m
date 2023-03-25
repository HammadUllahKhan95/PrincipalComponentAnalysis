clc;
clear;

%%
                 %%%%%%%% DATA DEFINTION %%%%%%%%
data = readtable('GEFCOMdata.csv',Delimiter=";");
data.Properties.VariableNames([1 2 3 4 5 6]) = {'Date','Hour','Zonal prices','System Load','Zonal Load','Day of Week'};

                        %%%%%%%% Calibration Window Parameters
mincol = zeros(size(data,1),1);
for i = 24:24:25968
    q = table2array(data(i-23:i,'Zonal prices'));
    s = min(q);
    t = ones(24,1);
    f = (s.*t);
    mincol(i+1:i+24) = f(1:24);
end
mincol = array2table(mincol);
mincol = mincol(25:25968+24,1);



datafinal  = [data mincol];
log_prices = log10(table2array(datafinal(:,"Zonal prices")));
log_min = log10(table2array(datafinal(1:25968,"mincol")));
log_load = log10(table2array(datafinal(:,"System Load")));

datafinal(:,"Zonal prices") = array2table(log_prices);
datafinal(:,"mincol") = array2table(log_min);
datafinal(:,"System Load") = array2table(log_load);

%%

da = table2array(data(:,"Date"));
cal = num2str(da);
dates= datetime(cal,'InputFormat','yyyyMMdd');
figure(1)
subplot(2,1,1);
plot(dates,table2array(data(:,"Zonal prices")));
hold on
plot([dates(17472,1) dates(17472,1)],[0 400],'k--')
hold on
plot([dates(18768,1) dates(18768,1)],[0 400],'k--')
xtickformat("dd.MM.yyyy")
xticks([dates(1,1) dates(366*24,1) dates(17473,1) dates(18769,1) dates(25968,1)])
xlim([dates(1,1) dates(25968,1)])
set(gca,'xticklabel',[])
ylabel('LMP [USD/MWh]')

load_array = table2array(data(:,"System Load"));
subplot(2,1,2);
plot(dates,load_array./1000);
hold on
plot([dates(17472,1) dates(17473,1)],[10 35],'k--')
hold on
plot([dates(18768,1) dates(18769,1)],[10 35],'k--')
xtickformat("dd.MM.yyyy")
xticks([dates(1,1) dates(366*24,1) dates(17473,1) dates(18769,1) dates(25968,1)])
xlim([dates(1,1) dates(25968,1)])
ylim([10 35])
set(gca,'XTickLabelRotation',45)
ylabel('System Load [GWh]')
%%
forecasts_all = zeros(8496,701);
for win = 1:701
    index_for_roll_cal = 0;
    T = 27+win;
    ind_d_entries = 0;
    vec = zeros(8496,1);
    index_for_data_entry = 1;
    datafinal = datafinal((1082-354-T)*24+1:end,:);
    for day = 1:354
        real_train = datafinal(index_for_roll_cal*24+1:(T+index_for_roll_cal)*24,:);
        for hour=1:24
            p = table2array(real_train(hour:24:end,"Zonal prices"));
            x = table2array(real_train(hour:24:end,"System Load"));
            k = table2array(real_train(hour:24:end,'mincol'));
            Dsat = table2array(real_train(hour:24:end, "Day of Week")) == 6;
            Dsun = table2array(real_train(hour:24:end, "Day of Week")) == 7;
            Dmon = table2array(real_train(hour:24:end, "Day of Week")) == 1;
            xx = table2array(datafinal(hour:24:end,"System Load"));
            DDsat = table2array(datafinal(hour:24:end, "Day of Week")) == 6;
            DDsun = table2array(datafinal(hour:24:end, "Day of Week")) == 7;
            DDmon = table2array(datafinal(hour:24:end, "Day of Week")) == 1;        
            % AR(1) - estimation
            pcal = p(1:T);
            xcal = x(1:T);
            kcal = k(1:T);
            Dsatcal = Dsat(1:T);
            Dsuncal = Dsun(1:T);
            Dmoncal = Dmon(1:T);
            yr= pcal(8:end); % for day d, d-1, d-2, d-7, load, min24hrbefore, 24th hour from day before, dummies  ...
            Xr = [ones(T-7,1) pcal(7:end-1) pcal(6:end-2) pcal(1:end-7) xcal(8:end) kcal(7:end-1) Dsatcal(8:end) Dsuncal(8:end) Dmoncal(8:end)];
            X_futr = [1 p(T) p(T-1) p(T-6) xx(T+1+ind_d_entries) k(T) DDsat(T+1+ind_d_entries) DDsun(T+1+ind_d_entries) DDmon(T+1+ind_d_entries)];
            % Regression, i.e., estimate betas
            beta = regress(yr,Xr);
            % Make prediction       
            pf1r = X_futr*beta;
            vec(index_for_data_entry) = pf1r; 
            index_for_data_entry = index_for_data_entry + 1;
        end
        index_for_roll_cal = index_for_roll_cal + 1;
        ind_d_entries = ind_d_entries + 1;
    end
    forecasts_all(:,win) = 10.^vec;
    datafinal  = [data mincol];
    log_prices = log10(table2array(datafinal(:,"Zonal prices")));
    log_min = log10(table2array(datafinal(1:25968,"mincol")));
    log_load = log10(table2array(datafinal(:,"System Load"))); 
    datafinal(:,"Zonal prices") = array2table(log_prices);
    datafinal(:,"mincol") = array2table(log_min);
    datafinal(:,"System Load") = array2table(log_load);
end
%%
writetable(array2table(forecasts_all),'thesisARX_results.csv');
%%
forecasts_all = readtable('thesisARX_results.csv');
forecasts_all = table2array(forecasts_all);
datafinal  = [data mincol];
%%
mae = [];
for i = 1:701
    mae(i,1) = mean(abs(forecasts_all(1297:end,i)-table2array(datafinal(18769:end,"Zonal prices"))));
end
plot(mae,'-')
xlim([0 701]);
ylabel('MAE')
xlabel('Calibration window length in days(T)')

%% 
%%%%%% WIN(28) %%%%%%%
mae_win_28 = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-forecasts_all(1297:end,1)));
%% 
%%%%%% WIN(364) %%%%%%%
mae_win_364 = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-forecasts_all(1297:end,337)));
%% 
%%%%%% AW(364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment364_728 = forecasts_all(:,[337,701]);
forecast_average364_728 = mean(forecasts_segment364_728,2);
mae_win_364_728 = mean(abs(data_need(1297:end) - forecast_average364_728(1297:end,1)));
%%
%%%%%% WAW(364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [337,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment364_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_364_728 = mean(abs(data_need(25+1296:end) - forecase_mae_final(1297:end,1)));
%%
%%%%%% expanding WAW(364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [337,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment364_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_364_728_exp =mean(abs(data_need(25+1296:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% WIN(728) %%%%%%%
mae_win_728 = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-forecasts_all(1297:end,701)));

%%
%%%%% AW(28:728) %%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_728 = forecasts_all(:,1:701);
forecast_average28_728 = mean(forecasts_segment28_728,2);
mae28_728 = mean(abs(data_need(1297:end) - forecast_average28_728(1297:end,1)));

%%
%%%%%% WAW(28:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:701
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,FC);
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:701) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_to_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:701
        data_forecasts = forecasts_all(1:24+indes*24,FC);
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:701) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_to_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%% AW(28:7:728) %%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_inc7_728 = forecasts_all(:,[1:7:701]);
forecast_average28_inc7_728 = mean(forecasts_segment28_inc7_728,2);
mae28_inc7_728 = mean(abs(data_need(1297:end) - forecast_average28_inc7_728(1297:end,1)));

%%
%%%%%% WAW(28:7:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:7:701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:101
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:101) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc7_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc7_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));;

%%
%%%%%% expanding WAW(28:7:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:7:701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:101
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:101) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc7_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc7_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%% AW(28:14:728) %%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_inc14_728 = forecasts_all(:,[1:14:701]);
forecast_average28_inc14_728 = mean(forecasts_segment28_inc14_728,2);
mae28_inc14_728 = mean(abs(data_need(1297:end) - forecast_average28_inc14_728(1297:end,1)));
%%
%%%%%% WAW(28:14:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:14:701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:51
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:51) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc14_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc14_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));
%%
%%%%%% expanding WAW(28:14:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:14:701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:51
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:51) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc14_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc14_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%% AW(28:28:728) %%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_inc28_728 = forecasts_all(:,[1:28:701]);
forecast_average28_inc28_728 = mean(forecasts_segment28_inc28_728,2);
mae28_inc28_728 = mean(abs(data_need(1297:end) - forecast_average28_inc28_728(1297:end)));
%%
%%%%%% WAW(28:28:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:28:701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:26
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:26) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc28_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc28_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28:28:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:28:701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:26
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:26) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc28_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc28_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% AW(56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment56_728 = forecasts_all(:,[29,701]);
forecast_average56_728 = mean(forecasts_segment56_728,2);
mae_win_56_728 = mean(abs(data_need(1297:end) - forecast_average56_728(1297:end)));
%%
%%%%%% WAW(56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [29,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment56_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_56_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));
%%
%%%%%% expanding WAW(56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [29,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment56_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_56_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% AW(28,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_and_728 = forecasts_all(:,[1,701]);
forecast_average28_and_728 = mean(forecasts_segment28_and_728,2);
mae_win_28_and_728 = mean(abs(data_need(1297:end) - forecast_average28_and_728(1297:end)));
%%
%%%%%% WAW(28,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_and_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_and_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:2
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:2) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_and_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_and_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%% AW(28:28:84,714:7:728) %%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_inc28_inc7 = forecasts_all(:,[1:28:57,687:7:701]);
forecast_average28_inc28_inc7 = mean(forecasts_segment28_inc28_inc7,2);
mae_inc28_inc7 = mean(abs(data_need(1297:end) - forecast_average28_inc28_inc7(1297:end)));
%%
%%%%%% WAW(28:28:84,714:7:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:28:57,687:7:701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:6
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:6) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc28_inc7(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc28_inc7 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28:28:84,714:7:728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1:28:57,687:7:701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:6
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:6) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_inc28_inc7(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_inc28_inc7_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% AW(28,56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_56_728 = forecasts_all(:,[1,29,701]);
forecast_average28_56_728 = mean(forecasts_segment28_56_728,2);
mae_win_28_56_728 = mean(abs(data_need(1297:end) - forecast_average28_56_728(1297:end)));
%%
%%%%%% WAW(28,56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:3
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:3) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28,56,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:3
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:3) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));


%%
%%%%%% AW(28,56,364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_56_364_728 = forecasts_all(:,[1,29,337,701]);
forecast_average28_56_364_728 = mean(forecasts_segment28_56_364_728,2);
mae_win_28_56_364_728 = mean(abs(data_need(1297:end) - forecast_average28_56_364_728(1297:end)));
%%
%%%%%% WAW(28,56,364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,337,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:4
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:4) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_364_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_364_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28,56,364,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,337,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:4
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:4) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_364_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_364_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% AW(28,56,721,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
forecasts_segment28_56_721_728 = forecasts_all(:,[1,29,694,701]);
forecast_average28_56_721_728 = mean(forecasts_segment28_56_721_728,2);
mae_win_28_56_721_728 = mean(abs(data_need(1297:end) - forecast_average28_56_721_728(1297:end)));
%%
%%%%%% WAW(28,56,721,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,694,701];
for d = 1:354
    data_actual = data_need(1+indes*24:24+indes*24,:);
    for FC = 1:4
        data_forecasts = forecasts_all(1+indes*24:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:4) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_721_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_721_728 = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%%%% expanding WAW(28,56,721,728)%%%%%%%
data_need  = table2array(datafinal(17473:end,"Zonal prices"));
iimae = [];
indes = 0;
weights_denom_mae=[];
num = [1,29,694,701];
for d = 1:354
    data_actual = data_need(1:24+indes*24,:);
    for FC = 1:4
        data_forecasts = forecasts_all(1:24+indes*24,num(1,FC));
        imae = 1/(mean(abs(data_actual - data_forecasts)));
        iimae(d,FC)= imae;
    end
    indes = indes + 1;
    weights_denom_mae(d,1) = sum(iimae(d,:));
end

weights_mae = [];
for dayy = 1:354
    weights_mae(dayy,1:4) = iimae(dayy,:)./weights_denom_mae(dayy,1);
end

weights_final_mae = repelem(weights_mae,24,1);
forecasts_mae_w_mult =forecasts_segment28_56_721_728(25:end,:).*weights_final_mae(1:end-24,:);
forecase_mae_final = sum(forecasts_mae_w_mult,2);
mae_inv__WAW_28_56_721_728_exp = mean(abs(data_need(1296+25:end) - forecase_mae_final(1297:end,1)));

%%
%%%PCA1
pca_tab = zeros(1320,701);
out_of_out_sample = forecasts_all(1297:end,:);
indexx = 0;
yy = zeros(1320,1);
mu_sig  = [];
finale = zeros(7200,1);
for j = 1:300
    part_of_forecast = out_of_out_sample(1+indexx*24:24+indexx*24,:);
    data_zzz = forecasts_all(1+indexx*24:1296+indexx*24,:);
    pca_tab(1:1296,:) = data_zzz;
    pca_tab(1297:1297+23,:) = part_of_forecast;
    data_actual_two = data_actual(1+indexx*24:1296+indexx*24,:);
    data_actual_two(1297:1320,1) = 0; 
    for i=1:1320
        mu = mean(pca_tab(i,:),2);
        sig = std(pca_tab(i,:),0);
        if i >1296
            mu_sig(i-1296,1:2) = [mu sig];
        end
        pca_tab_f(i,:) = (pca_tab(i,:)-mu)./sig;
        yy(i,1) = (data_actual_two(i,:)-mu)./sig;
    end
    yy = yy(1:1296,1);
    [L_nn,~] = eigs(pca_tab_f'*pca_tab_f/54,1);
    L_nn = sqrt(701)*L_nn;
    F_tt = 1/((54+1)*24)*pca_tab_f*L_nn;
    n = length(54*24);
    betas = regress(yy(1:1296,1),F_tt(1:1296,:));
    X_future = F_tt(1297:end,:);
    vals = X_future*betas;
    final = (vals.*mu_sig(:,2))+mu_sig(:,1);
    finale(1+indexx*24:24+indexx*24,1) = final; 
    indexx = indexx + 1;
end
final_pca1_mae = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-finale));




%%
%%%PCA2
pca_tab = zeros(1320,701);
out_of_out_sample = forecasts_all(1297:end,:);
indexx = 0;
yy = zeros(1320,1);
mu_sig  = [];
finale = zeros(7200,1);
for j = 1:300
    part_of_forecast = out_of_out_sample(1+indexx*24:24+indexx*24,:);
    data_zzz = forecasts_all(1+indexx*24:1296+indexx*24,:);
    pca_tab(1:1296,:) = data_zzz;
    pca_tab(1297:1297+23,:) = part_of_forecast;
    data_actual_two = data_actual(1+indexx*24:1296+indexx*24,:);
    data_actual_two(1297:1320,1) = 0; 
    for i=1:1320
        mu = mean(pca_tab(i,:),2);
        sig = std(pca_tab(i,:),0);
        if i >1296
            mu_sig(i-1296,1:2) = [mu sig];
        end
        pca_tab_f(i,:) = (pca_tab(i,:)-mu)./sig;
        yy(i,1) = (data_actual_two(i,:)-mu)./sig;
    end
    yy = yy(1:1296,1);
    [L_nn,~] = eigs(pca_tab_f'*pca_tab_f/54,2);
    L_nn = sqrt(701)*L_nn;
    F_tt = 1/((54+1)*24)*pca_tab_f*L_nn;
    n = length(54*24);
    betas = regress(yy(1:1296,1),F_tt(1:1296,:));
    X_future = F_tt(1297:end,:);
    vals = X_future*betas;
    final = (vals.*mu_sig(:,2))+mu_sig(:,1);
    finale(1+indexx*24:24+indexx*24,1) = final; 
    indexx = indexx + 1;
end
final_pca2_mae = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-finale));

%%

%%%PCA3
pca_tab = zeros(1320,701);
out_of_out_sample = forecasts_all(1297:end,:);
indexx = 0;
yy = zeros(1320,1);
mu_sig  = [];
finale = zeros(7200,1);
for j = 1:300
    part_of_forecast = out_of_out_sample(1+indexx*24:24+indexx*24,:);
    data_zzz = forecasts_all(1+indexx*24:1296+indexx*24,:);
    pca_tab(1:1296,:) = data_zzz;
    pca_tab(1297:1297+23,:) = part_of_forecast;
    data_actual_two = data_actual(1+indexx*24:1296+indexx*24,:);
    data_actual_two(1297:1320,1) = 0; 
    for i=1:1320
        mu = mean(pca_tab(i,:),2);
        sig = std(pca_tab(i,:),0);
        if i >1296
            mu_sig(i-1296,1:2) = [mu sig];
        end
        pca_tab_f(i,:) = (pca_tab(i,:)-mu)./sig;
        yy(i,1) = (data_actual_two(i,:)-mu)./sig;
    end
    yy = yy(1:1296,1);
    [L_nn,~] = eigs(pca_tab_f'*pca_tab_f/54,3);
    L_nn = sqrt(701)*L_nn;
    F_tt = 1/((54+1)*24)*pca_tab_f*L_nn;
    n = length(54*24);
    betas = regress(yy(1:1296,1),F_tt(1:1296,:));
    X_future = F_tt(1297:end,:);
    vals = X_future*betas;
    final = (vals.*mu_sig(:,2))+mu_sig(:,1);
    finale(1+indexx*24:24+indexx*24,1) = final; 
    indexx = indexx + 1;
end
final_pca3_mae = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-finale));

%%

%%%PCA4
pca_tab = zeros(1320,701);
out_of_out_sample = forecasts_all(1297:end,:);
indexx = 0;
yy = zeros(1320,1);
mu_sig  = [];
finale = zeros(7200,1);
for j = 1:300
    part_of_forecast = out_of_out_sample(1+indexx*24:24+indexx*24,:);
    data_zzz = forecasts_all(1+indexx*24:1296+indexx*24,:);
    pca_tab(1:1296,:) = data_zzz;
    pca_tab(1297:1297+23,:) = part_of_forecast;
    data_actual_two = data_actual(1+indexx*24:1296+indexx*24,:);
    data_actual_two(1297:1320,1) = 0; 
    for i=1:1320
        mu = mean(pca_tab(i,:),2);
        sig = std(pca_tab(i,:),0);
        if i >1296
            mu_sig(i-1296,1:2) = [mu sig];
        end
        pca_tab_f(i,:) = (pca_tab(i,:)-mu)./sig;
        yy(i,1) = (data_actual_two(i,:)-mu)./sig;
    end
    yy = yy(1:1296,1);
    [L_nn,~] = eigs(pca_tab_f'*pca_tab_f/54,4);
    L_nn = sqrt(701)*L_nn;
    F_tt = 1/((54+1)*24)*pca_tab_f*L_nn;
    n = length(54*24);
    betas = regress(yy(1:1296,1),F_tt(1:1296,:));
    X_future = F_tt(1297:end,:);
    vals = X_future*betas;
    final = (vals.*mu_sig(:,2))+mu_sig(:,1);
    finale(1+indexx*24:24+indexx*24,1) = final; 
    indexx = indexx + 1;
end
final_pca4_mae = mean(abs(table2array(datafinal(18769:end,"Zonal prices"))-finale));
