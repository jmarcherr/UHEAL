function [age,aud_L,aud_R,aud_freq,gender] = get_aud(dataalm)
%Get age and audiograms for participants
    aud = dataalm.aud;
    aud_freq = squeeze(aud(:,1,1));
    aud_L = squeeze(aud(:,2,1));
    aud_R = squeeze(aud(:,2,2));
    % Age
    str=dataalm.per.BirthDate.Text;
    AUDdat = dataalm.per.AUDdate(1:10);

    birth_numdate=datenum(str,'YYYY-mm-DD');
    %AUD_date = datenum(AUDdat,'YYYY-mm-DD');
    if isempty(dataalm.rds)
        age =dataalm.subinfo.age;
    else
    AUD_date = datenum(dataalm.rds.date);
    age=datestr(AUD_date-birth_numdate,'YYYY-mm-DD');
    age = str2num(age(1:4));
    end
        if strcmp(dataalm.id,'UH091')
            age=19;
        end
    if strcmp(dataalm.per.Gender.Text,'Female')
        gender = 1;
    else
        gender =2;
    end
end

