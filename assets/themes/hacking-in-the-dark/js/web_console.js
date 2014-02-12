APP_TITLE = 'WebConsole'
VERSION = '0.1.5';
AUTHOR = 'Jonathan Tsai';
AUTHOR_EMAIL = 'akajontsai-devel@yahoo.com';
COPYRIGHT_YEAR = 2013;

YUI.add('web-console', function (Y) {
    // CSS selectors
    var CSS_CLASS_CONSOLE = 'console'

    var CSS_ID_WEBCONSOLE_PANEL = 'webconsole_panel';

    // "Local" globals
    KEY_MAP = Y.Node.DOM_EVENTS.key.eventDef.KEY_MAP;

    var HTML_NEWLINE = '\n';
    var PROMPT = '>>> ';

    // member variables
    var _instance = null;

    // member functions
    function WebConsole(config) {
        var _panel;
        var commandHistory = [];

        this.test = function() {
            Y.log('WebConsole.test()');
        };

        this.keyboardShortcutHandlers = {
            'enter' : function(e) {
                e.preventDefault();
                var input = _instance.getConsoleInput();
                _instance.processCommand(input);
                _instance.cprompt();
                _instance.cend();
            },
            'backspace' : function(e) {
                _instance.cend();
                var input = _instance.getConsoleInput();
                if (input === '') {
                    // don't delete any characters if at beginning of line
                    e.preventDefault();
                }
            }
        };

        this.consoleCommands = {
            'help' : function(args) {
                if (args.length === 0) {
                    _instance.cout('Type "help commands" to see a list of commands');
                    _instance.cline();
                    _instance.cout('What would you like help with?');
                } else {
                    if (args[0] === 'commands') {
                        var commands = [];
                        for (var command in _instance.consoleCommands) {
                            commands.push(command);
                        }
                        commands.sort();
                        _instance.cout('Commands: ' + commands.join(', '));
                    } else {
                        _instance.cout('No help for: ' + args.join(' '));
                    }
                }
            },
            'copyright' : function(args) {
                _instance.cout(Y.Lang.sub('Copyright {copyright_year} by {author} <{author_email}>', { copyright_year : COPYRIGHT_YEAR, author : AUTHOR, author_email : AUTHOR_EMAIL }));
            },
            'credits' : function(args) {
                _instance.cout(Y.Lang.sub('Created by {author}. Inspired by Python, Unix, MUDs, and countless other terminals.', { author : AUTHOR }));
            },
            'license' : function(args) {
                _instance.cout('Type license() to see the full license text');
            },
            'license()' : function(args) {
                _instance.cout('MIT licensed: https://github.com/jontsai/jekyll-theme-hacking-in-the-dark/blob/master/LICENSE');
            },
            'exit' : function(args) {
                _instance.hide();
            },
            'quit' : function(args) {
                _instance.consoleCommands['exit'](args);
            },
            'clear' : function(args) {
                _instance.clear();
            },
            'cls' : function(args) {
                _instance.consoleCommands['clear'](args);
            },
            'history' : function(args) {
                for (var i=0; i < commandHistory.length; ++i) {
                    var command = commandHistory[i];
                    _instance.cout(i + ': ' + command);
                }
            },
            'echo' : function(args) {
                _instance.cout(args.join(' '));
            },
            'print' : function(args) {
                _instance.consoleCommands['echo'](args);
            },
            'time' : function(args) {
                var d = new Date();
                //_instance.cout(d.toISOString());
                _instance.cout(d.toLocaleString());
            },
            'utc' : function(args) {
                var d = new Date();
                _instance.cout(d.toUTCString());
            },
            'agent' : function(args) {
                var userAgent = navigator.userAgent;
                var appVersion = navigator.appVersion;
                _instance.cout(userAgent + ' ' + appVersion);
            },
            'referrer' : function(args) {
                var referrer = document.referrer;
                _instance.cout(referrer);
            },
            'ip' : function(args) {
                var ip;
                if (typeof(java) !== 'undefined' && typeof(java.net) !== 'undefined') {
                    ip = '' + java.net.InetAddress.getLocalHost().getHostAddress();
                } else {
                    ip = 'unknown';
                }
                _instance.cout('Your IP address as detected by JavaScript: ' + ip);
            },
            'location' : function(args) {
                if (navigator.geolocation) {
                    _instance.cout('Getting GPS location...');
                    navigator.geolocation.getCurrentPosition(function(position) {
                        var lat = position.coords.latitude;
                        var lng = position.coords.longitude;
                        _instance.cout('GPS: ' + lat + ', ' + lng);
                        _instance.cprompt();
                        if (confirm('Open in Google Maps?')) {
                            window.open('http://maps.google.com/maps?q=' + lat + ',' + lng);
                        }
                    });
                } else {
                    _instance.cout('navigator.geolocation is not available');
                }
            },
            'curl' : function(args) {
                if (args.length !== 1) {
                    _instance.cout('invalid URI: ' + args.join(' '));
                } else {
                    var uri = args[0];
                    var cfg = {
                        on : {
                            complete : function(transactionId, response, args) {
                                var responseText = response.responseText;
                                _instance.cout(responseText);
                            },
                            failure : function() {
                            }
                        }
                    }
                    var response = Y.io(uri, cfg);
                }
            },
            '!' : function(args) {
                var previousCommand = commandHistory[commandHistory.length - 2];
                _instance.processCommand(previousCommand);
            }
        };

        /**
         * getPanel
         *
         * Gets the Panel containing the web console
         * Initializes a new Panel if it does not exist yet
         */
        this.getPanel = function() {
            if (typeof _panel === 'undefined') {
                var panelContent = Y.Node.create(Y.one('#' + CSS_ID_WEBCONSOLE_PANEL).getHTML());
                var documentBody = Y.one(document.body);
                var winWidth = documentBody.get('winWidth');
                var winHeight = documentBody.get('winHeight');

                var panelCfg = {
                    srcNode : panelContent,
                    width : winWidth * 0.8,
                    height : winHeight * 0.8,
                    zIndex : 90,
                    modal : true,
                    visible : false,
                    render : true,
                    hideOn : []
                };

                var panel = new Y.Panel(panelCfg);

                _panel = panel;
                // set up event delegation for key events in panel
                var boundingBox = panel.get('boundingBox');
                boundingBox.delegate('clickoutside', function(e) { _instance.hide(); }, '*');
                boundingBox.delegate('key', function(e) { _instance.hide(); }, 'down:esc', '*');

                var headerContentParams = {
                    app_title: APP_TITLE,
                    version: VERSION,
                    author: AUTHOR
                }
                var headerContent = Y.Lang.sub('{app_title} {version} by {author}', headerContentParams);
                var header = boundingBox.one('.yui3-widget-hd');
                header.setHTML(headerContent);

                this.initConsole();
            }
            return _panel;
        };

        /**
         * initConsole
         *
         * initializes the WebConsole text, etc
         */
        this.initConsole = function() {
            var introParams = {
                app_title: APP_TITLE,
                version: VERSION,
                author: AUTHOR,
                author_email: AUTHOR_EMAIL,
                copyright_year: COPYRIGHT_YEAR
            };
            var introMessage = Y.Lang.sub('{app_title} {version} by {author} <{author_email}> (c) {copyright_year}', introParams);
            this.cout(introMessage, false);
            this.cout('Type "help", "copyright", "credits" or "license" for more information.');
            this.cprompt();
            this.cend();

            var console = this.getConsole();
            for (var keyName in this.keyboardShortcutHandlers) {
                var handler = this.keyboardShortcutHandlers[keyName];
                console.delegate('key', handler, 'down:' + keyName, '*');
            }
        };

        /**
         * getConsole
         *
         * @returns the console DOM node
         */
        this.getConsole = function() {
            var console = _panel.get('boundingBox').one('.' + CSS_CLASS_CONSOLE);
            return console;
        };

        /**
         * getConsoleInput
         *
         * @returns the last line of input from the console
         */
        this.getConsoleInput = function() {
            var console = this.getConsole();
            var lines = console.get('value').split('\n');
            var lastLine = lines[lines.length - 1];
            var input;
            if (lastLine.indexOf(PROMPT) === 0) {
                input = Y.Lang.trim(lastLine.split(PROMPT)[1]);
            } else {
                input = Y.Lang.trim(lastLine);
            }
            return input;
        };

        /**
         * processCommand
         *
         */
        this.processCommand = function(input) {
            commandHistory.push(input);

            var tokens = input.split(' ');
            var command = tokens.shift();
            if (command !== '') {
                var args = tokens.length > 0? tokens : [];
                var handler = this.consoleCommands[command];
                if (typeof handler === 'function') {
                    handler(args);
                } else {
                    _instance.cout(command + ': command not found');
                }
            }
        };

        /**
         * cout
         *
         * Prints the given string to console
         */
        this.cout = function(s, addNewline) {
            addNewline = addNewline === undefined ? true : addNewline;
            var console = this.getConsole();
            var value = console.get('value');
            if (addNewline) {
                value += HTML_NEWLINE;

            }
            value += s;
            console.set('value', value);
        };

        /**
         * cline
         * Prints an empty line to the console
         */
        this.cline = function() {
            this.cout('');
        }

        /**
         * cprompt
         *
         * Prints the given string to console, prepended with HTML_NEWLINE + PROMPT
         */
        this.cprompt = function(s) {
            s = s === undefined? '' : s;
            this.cout(HTML_NEWLINE + PROMPT + s, false);
        };

        /**
         * cend
         *
         * Focuses on the console and puts cursor at end of text
         */
        this.cend = function() {
            var console = this.getConsole();
            console.focus();
            var value = console.get('value');
            console.set('value', '');
            console.set('value', value);
        };

        /**
         * clear
         *
         * Clears the console
         */
        this.clear = function() {
            var console = this.getConsole();
            console.set('value', '');
        }

        /**
         * show
         *
         * shows the WebConsole panel
         */
        this.show = function() {
            var panel = this.getPanel();
            var boundingBox = panel.get('boundingBox');
            var width = boundingBox.getStyle('width');
            panel.show();
            boundingBox.setStyle('top', '30px');
            boundingBox.setStyle('left', '-' + width);
            boundingBox.transition({
                duration: 0.25,
                left: '0px'
            });
            this.cend();
        };

        /**
         * hide
         *
         * hides the WebConsole panel
         */
        this.hide = function() {
            var panel = this.getPanel();
            var boundingBox = panel.get('boundingBox');
            var width = boundingBox.getStyle('width');
            boundingBox.transition({
                duration: 0.25,
                left: '-' + width
            }, function() {
                panel.hide();
            });
        };

        /**
         * toggle
         *
         * toggles the WebConsole panel shown/hidden state
         */
        this.toggle = function() {
            var panel = this.getPanel();
            if (panel.get('visible')) {
                this.hide();
            } else {
                this.show();
            }
        }
    }

    function getInstance() {
        if (!_instance) {
            var cfg = {};
            _instance = new WebConsole(cfg);
        }
        var console = _instance;
        return console;
    }

    function toggle() {
        var console = getInstance();
        console.toggle();
    }

    function isExpanded() {
        var console = getInstance();
        var expanded = true && console.getPanel().get('visible');
        return expanded;
    }

    Y.WebConsole = {
        toggle : toggle,
        isExpanded : isExpanded
    }
}, VERSION, {
    requires: [
        'node',
        'event',
        'event-key',
        'transition',
        'io'
    ]
});
