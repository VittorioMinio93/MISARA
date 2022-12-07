%%%%%%%%%%%%% Get and set the section of the menu panel%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio

function menu_cbk(~,eventdata)

old = get(eventdata.OldValue,'Userdata');%Get the old section value of the menu panel
new = get(eventdata.NewValue,'Userdata');%Get the new section value of the menu panel

set (old, 'visible', 'off')%Set the old section value of the menu panel as off
set (new, 'visible', 'on')%Set the new section value of the menu panel as on

end