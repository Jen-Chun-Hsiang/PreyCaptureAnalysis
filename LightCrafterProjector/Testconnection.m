
host = '10.26.29.138';
port = 443;

% host = '127.0.0.1';
% port = 8080;
%%
version('-java')
socket = java.net.Socket();
addr = java.net.InetSocketAddress(host, port);
socket.connect(addr, 135)

%%

connection = netbox.Connection(host, port);

%%
port = 12345;
server = tcpserver(host, port);
client = tcpclient(host, port);
data = 'hello';
write(client, data, 'string');

if server.NumBytesAvailable >0 
    received_data = read(server, server.NumBytesAvailable, 'string');
    disp(['server receive:' received_data]);
end

%%
stageClient = stage.core.network.StageClient();
stageClient.connect(host, port);
stageClient.getCanvasSize
stageClient.disconnect()