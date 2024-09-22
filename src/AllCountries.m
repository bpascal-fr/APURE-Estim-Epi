% List all the 200+ countries for which  Covid19 new infections counts 
% are available in the Johns Hopkins University repository: 
% https://coronavirus.jhu.edu/.
%
% B. Pascal,
% March, 2024


function AllCountries(display)

    % Inputs:  - display: 1 for displaying the list of all available countries,
    %                     0 does noting.

    if display

        % Load data made available by Johns Hopkins University from folder data
        outfilename = 'data/COVID-19-JHU_stored.csv';

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
            datesAscii              = T.Properties.VariableNames ;
            Dates                   = datenum(datesAscii(5:end));
        end


        % Store data in correct format for processing
        CountriesTMP   = table2array(T(:,2)) ;
        Provinces      = table2array(T(:,1)) ;
        TConfirmedTMP  = table2array(T(:,5:end));
        Compteur       = 1 ;
        k              = 2 ;

        while k <= length(CountriesTMP)
            All_Countries{Compteur}            = CountriesTMP{k} ;
            if strcmp(Provinces{k},'')
                TConfirmed(Compteur,:)     = TConfirmedTMP(k,:) ;
                k                          = k+1 ;
            else
                TConfirmed(Compteur,:)     = TConfirmedTMP(k,:) ;
                k2                         = 1 ;
                while strcmp(CountriesTMP{k+k2},CountriesTMP{k})
                    TConfirmed(Compteur,:) = TConfirmed(Compteur,:) + TConfirmedTMP(k+k2,:) ;
                    k2                     = k2+1 ;
                end
                k                          = k+k2 ;
            end
            Compteur                       = Compteur + 1 ;
        end

        disp(cell2table(All_Countries','VariableNames',"Country name"))

    end

    

end