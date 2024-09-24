% Load Covid19 daily new infections counts in 200+ countries worldwide during
% ~ 3 years of Covid19 pandemic made available by the Johns Hopkins
% University on the repository: https://coronavirus.jhu.edu/ from confirmed
% cases collected by National Health Agencies.
%
% Then turn the daily infection counts into weekly infection counts by 
% aggregating counts over seven days in each territory separately and 
% compute the associated  weekly global infectiousness.
% 
% The daily discretized serial interval function is also weekly discretized
% using the left-rectangle integration method with constant step size of 
% one day length.
%
%
% Last release of JHU repository on October 3, 2023.
%
% On March 31, 2024, data are still available at https://coronavirus.jhu.edu/
% the present code proposes to either:
% - download them from the Internet,
% - load a .csv file stored in the folder data.
%
% See data/README.md for more detail.
%
% References
%
% - Nash, R. K., Bhatt, S., Cori, A., & Nouvellet, P. (2023). Estimating the 
% epidemic reproduction number from temporally aggregated incidence data: 
% A statistical modelling approach and software tool. 
% PLOS Computational Biology, 19(8), e1011439.
%
% - Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a
% Nonstationary Autoregressive Model for Data-Driven Reproduction Number
% Estimation. Preprint. arXiv:2409.14937.
%
% B. Pascal and P. Abry, March 2024.



function [Z_Week, Phi_Z_Week, M_Week] = load_JHU_Weekly(User_Country,opts)

    % Inputs:  - User_Country: strings containing the names of the country which is to be monitored, 
    %                 e.g., User_Country = "France".
    %          - opts: structure indicating the time period selected by the user containing (optional)
    %                   opts.LastDay: last day of the time period in format 'YYYY-MM-DD' (default '2023-03-09', latest day possible)
    %                   opts.W: number of weeks of the time period (default: total number of weeks available in JHU repository).
    %                   opts.Download: 0 for loading the .csv in folder data, 1 for downloading data from https://coronavirus.jhu.edu/
    %                   opts.Phi: daily discretized serial interval function to be used (default: daily discretized Covid19 serial interval function)
    %                   opts.FontSize: desired FontSize for the plots (optional)
    %
    %
    % Outputs: - Z_Week: weekly aggregated infection counts
    %          - Phi_Z_Week: associated global infectiousness
    %          - M_Week: weekly and daily data and model parameters, structure containing
    %                   M_Week.Phi: weekly discretized serial interval function
    %                   M_Week.Phi_Day: daily discretized serial interval function
    %                   M_Week.Dates: dates corresponding to the last days of the weeks over which infection counts are aggreagted.
    %                   M_Week.Dates_Day: dates corresponding to all days of the considered time period.
    %                   M_Week.Z_day: daily new infection counts during the whole time period.
    %                   M_Week.Phi_Z_day: daily infectiousness during the whole time period.
    %                   M_Week.Y_SY: weekly new infection counts during the extended time period.



    if nargin < 2

        % Last day of the time period
        ChosenTime = [2023,03,09];

        % Length of the studied period
        TUser      = -1;

        % Load data from local folder
        Download   = 0;

        % Fontsize of the plot
        FontSize   = 22.5 ;

    else
        
        if ~isfield(opts,'LastDay'),  opts.LastDay  = '2023-03-09'; end
        if ~isfield(opts,'W'),        opts.W        = -2; end
        if ~isfield(opts,'Download'), opts.Download = 0; end
        if ~isfield(opts,'FontSize'), opts.FontSize = 22.5; end
        
        % Last day of the period
        LastDay    = opts.LastDay;
        ChosenTime = [str2double(LastDay(1:4)),str2double(LastDay(6:7)),str2double(LastDay(9:10))];

        % Length of the studied period in days
        TUser      = 7 * (opts.W + 1) + 1;

        % Load or download data
        Download   = opts.Download;

        % Fontsize of the plot
        FontSize   = opts.FontSize;

    end
   

    %% LOAD DATA MADE AVAILABLE BY JOHNS HOPKINS UNIVERSITY

    % Download or load from stored data
    if Download
        url         = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv';
        filename    ='COVID-19-JHU_downloaded.csv';
        outfilename = websave(filename,url);
    else
        outfilename = 'data/COVID-19-JHU_stored.csv';
    end


    % Handle different encoding of the dates depending on Matlab version
    [~,d_MAT]      = version;      % Matlab release date
    y_MAT          = d_MAT(end-3:end); % year of release


    % Load Covid19 confirmed cases in all countries
    if strcmp(y_MAT, 2019)
        T       = readtable(outfilename);
        Dates   = datenum(table2array(T(1,5:end)));
        T       = readtable(outfilename,'ReadVariableNames',true);
    else
        opts                    = detectImportOptions(outfilename);
        opts.VariableNamesLine  = 1;
        opts.VariableNamingRule = 'preserve';
        T                       = readtable(outfilename,opts);
        datesAscii              = T.Properties.VariableNames;
        Dates                   = datenum(datesAscii(5:end));
    end


    % Store data in correct format for processing
    CountriesTMP   = table2array(T(:,2));
    Provinces      = table2array(T(:,1));
    TConfirmedTMP  = table2array(T(:,5:end));
    Compteur       = 1;
    k              = 2;

    while k <= length(CountriesTMP)
        All_Countries{Compteur}            = CountriesTMP{k};
        if strcmp(Provinces{k},'')
            TConfirmed(Compteur,:)     = TConfirmedTMP(k,:);
            k                          = k+1;
        else
            TConfirmed(Compteur,:)     = TConfirmedTMP(k,:);
            k2                         = 1;
            while strcmp(CountriesTMP{k+k2},CountriesTMP{k})
                TConfirmed(Compteur,:) = TConfirmed(Compteur,:) + TConfirmedTMP(k+k2,:);
                k2                     = k2+1;
            end
            k                          = k+k2;
        end
        Compteur                       = Compteur + 1;
    end

    %% STORE INDICES CORRESPONDING TO COUNTRIES REQUIRED BY THE USER
    % If the country is not found in the list or misspelled it is ignored and
    % not included in the list of analyzed countries

    N_Countries = [];
    Countries   = [];

    for n = 1:length(User_Country)

        N_n            = find(strcmp(All_Countries,User_Country(n)));

        if isempty(N_n)

            error(strcat(User_Country(n),' not found. Look at the list of available countries provided by calling AllCountries(1).'))

        else

            N_Countries = [N_Countries, N_n];
            Countries   = [Countries, User_Country(n)];

        end


    end


    %% EXTRACT NEW INFECTION COUNT MULTIVARIATE TIME SERIES

    % Daily cumulated number of cases
    Confirmed = TConfirmed(N_Countries,:);

    % Daily new infection counts
    Zdata     = diff(Confirmed,1,2);
    Dates     = Dates(2:end);

    % Store dates
    aa        = datevec(Dates(end)); aa = aa(1:3);
    bb        = datevec(Dates(1))  ; bb = bb(1:3);
    tdate     = datetime(bb(1),bb(2),bb(3)):datetime(aa(1),aa(2),aa(3));

    % Start at the beginning of the epidemic for the studied countries
    begin_countries = zeros(size(Countries));

    for n = 1:length(Countries)

        ind_begin          = find(Zdata(n,:) > 0);
        begin_countries(n) = ind_begin(1);

    end

    begin_pandemic  = min(begin_countries);

    Zdata           = Zdata(:,begin_pandemic:end);
    tdate           = tdate(begin_pandemic:end);
    Dates           = Dates(begin_pandemic:end);

    %% PREPROCESSING

    % Remove negative counts
    Zdata(Zdata < 0) = 0;

    % Apply a sliding median to discard outlier values
    alpha            = 7;
    Zalpha           = daily_sliding_median(Zdata,alpha,7);

    %% COMPUTE GLOBAL INFECTIOUSNESS

    % Serial interval function modeled by a Gamma distribution with
    % - mean: 6.6 days
    % - standard deviation 3.5 days
    % truncated at 25 days and normalized.
    shape    = 1/0.28;
    scale    = 1.87;
    tau_phi  = 25;
    Phi      = gampdf(0:tau_phi,shape,scale);
    Phi      = Phi/sum(Phi); % normalize the weights applied to past tau_phi infection counts

    % Infectiousness: weighted sum of past infection counts
    [Phi_Z,Zalpha] = Psi_normal(Zalpha,Phi);
    tdate          = tdate(2:end);
    Dates          = Dates(2:end);
    % one day cropped because infectiousness not defined at day 1 due to border effect induced by absence of any past count

    %% FOCUS ON THE SPECIFIED TIME PERIOD

    % Find the TUser dates
    tChosenTime = datenum(ChosenTime);
    index       = find(Dates == tChosenTime);

    % Define last day
    if isempty(index)
        warning('Required last day not covered by JHU repository, replaced by the last availble day: March 9, 2023.')       
        index = size(Zalpha,2);
    end

    % Define the length of the time period
    if TUser < 0
        date = 1:index;
    else
        if index-TUser+1 <= 0
            warning(strcat('Time period not fully included in JHU repository hence adjusted to T = '," ",num2str(index),' days.'))
            date        = 1:index;
        else
            date        = index-TUser+1:index;
        end
    end

    %% INFECTION COUNTS, INFECTIOUSNESS AND DATES

    % Extract infection counts, infectiousness and dates for the specified time period
    Z_Day          = Zalpha(:,date);
    Phi_Z_Day      = Phi_Z(:,date);
    Dates_Day      = tdate(date);

    %% TURN THE DAILY DATA INTO WEEKLY DATA

    % Aggregation of daily infection counts over seven days to get weekly counts
    W_User          = floor(length(Z_Day)/7);
    Z_Week          = zeros(size(Z_Day,1),W_User);
    for w = 1:W_User
        tmp_Z         = Z_Day(:,7*(w-1)+1:7*w);
        Z_Week(:,w)   = sum(tmp_Z,2);
    end

    % Corresponding dates
    Dates_Day    = Dates_Day(1:7*W_User); 
    Dates_Week   = Dates_Day(7:7:end);

    % Construct the discretized serial interval function
    if ~isfield(opts,'Phi') % default: Covid19 serial interval function
        tau_day  = 25;                            % memory horizon
        shape    = 1/0.28;                        % shape parameter
        scale    = 1.87;                          % scale parameter
        Phi      = gampdf(0:tau_day,shape,scale); % discretized Gamma probability density function
    else
        % the serial interval function should start with a zero and should be a row vector
        if ~(opts.Phi(1) == 0)
            tau_day = length(opts.Phi);
            Phi     = reshape(opts.Phi,1,tau_day);
            Phi     = [0, Phi];
        else
            tau_day = length(opts.Phi) - 1;
            Phi     = reshape(opts.Phi,1,tau_day+1);
        end
    end
    % tau - 1 should be a multiple of seven
    if mod(tau_day,7)
        add_day = 7 - mod(tau_day,7);
        tau_day = tau_day + add_day;
        Phi     = [Phi, zeros(1,add_day)];
    end
    Phi         = Phi/sum(Phi);             % normalize the serial interval function (optional)

    % Aggregation over seven days to get a weekly discretized interval function
    Phi_Week = zeros(1,tau_day/7);
    for w    = 1:tau_day/7
        tmp_Phi         = Phi(7*(w-1)+2:7*w+1);
        Phi_Week(w+1)   = sum(tmp_Phi);
    end

    % Store the infection counts during the extended time period
    M_Week.Y_SY         = Z_Week;

    % Weekly global infectiousness
    [Phi_Z_Week,Z_Week] = Psi_normal(Z_Week,Phi_Week);

    % One week cropped because infectiousness not defined at week 1 due to border effect induced by absence of any past count
    Z_Day               = Z_Day(:,8:end);
    Dates_Day           = Dates_Day(8:end);
    Dates_Week          = Dates_Week(2:end);

    %% STORE THE DAILY DATA AND PARAMETERS OF THE MODEL

    % Weekly and daily serial interval functions
    M_Week.Phi         = Phi_Week;
    M_Week.Phi_Day     = Phi;

    % Dates on a weekly or daily basis
    M_Week.Dates       = Dates_Week;
    M_Week.Dates_Day   = Dates_Day;

    % Daily infection counts and infectiousness
    M_Week.Z_Day       = Z_Day;
    M_Week.Phi_Z_Day   = Phi_Z_Day;

    %% DISPLAY 

    % Messages to the user
    monitored = User_Country(1); 
    for n     = 2:length(User_Country)
        monitored = strcat(monitored,", ",User_Country(n));
    end

    disp('---------------------------------------------------------')
    disp('COUNTRIES AND TIME PERIOD CHOSEN FOR ANALYSIS:')
    disp(strcat("Monitored countries: ",monitored,"."))
    disp(strcat("Time period: from ",string(Dates_Day(1))," to ",string(Dates_Day(end))))
    disp('---------------------------------------------------------')

    f2              = figure(2); clf
    p               = plot(Dates_Week, Z_Week,'-','linewidth',2,'color','black') ;
    grid on ; hold on
    q               = plot(Dates_Week, Phi_Z_Week,'-.','linewidth',2,'color','black') ;
    leg             = legend([p,q],'$\mathrm{Z}_t$','$\Phi_t(\mathbf{Z})$','location','best');
    leg.Interpreter = 'Latex';
    leg.Color       = 'none';
    leg.FontSize    = FontSize;
    VX              = [Dates_Day(1) Dates_Day(end)] ;
    xlim(VX)
    title(strcat("Weekly aggregated COVID-19 infection counts in ",User_Country),'Interpreter','Latex')
    set(gca,'ticklabelinterpreter','Latex','fontsize',FontSize,'color','None')
    f2.Position     = [141 329 1033 314];

end