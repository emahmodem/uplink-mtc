function [response] = send_report_via_email(subject,report)

% Sending Email
sending_email = '';
mail_server = '';
user_name = '';
password = '';


receiving_email = '';


% Set Sending Email credentials
setpref('Internet','E_mail',sending_email);
setpref('Internet','SMTP_Server',mail_server);
setpref('Internet','SMTP_Username',user_name);
setpref('Internet','SMTP_Password',password);

% Set security properties
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');


sendmail(receiving_email,subject,report) ;

response = 'done';
end
