# Einführung in zeitkontinuierliche Übertragungsfunktionen

## Überblick
Dies ist ein Vortrag über zeitkontinuierliche Übertragungsfunktionen. Als Teil dessen werden die Übertragungsfunktionen $G(s)$ und die Bode-Diagramme $G(j\omega)$ vorgestellt. Es werden einige typische Beispiele aus Lehrbüchern vorgestellt und die mathematische Formulierung der Übertragungsfunktionen erläutert.

## Installation

Die HTML-Dateien können direkt geöffnet werden. Laden Sie sich eine ZIP-Datei herunter und entpacken Sie sie. Ausgepackt finden Sie die HTML-Datei `vortrag.html` und in ausgepackter Form funktionieren die Verweise und Links auf die weiteren Unterseiten.

## Editieren des Vortrags

Wenn Sie den Vortrag editieren möchten, dann müssen Sie die Quell-Datei `vortrag.qmd` editieren. Dazu benötigen Sie Quarto und eine Python-Installation (die Präsentation enthält Grafiken, die über Python generiert werden). Weitere Informationen über das Präsentationsformat finden Sie in der [Quarto-Dokumentation](https://quarto.org/docs/guide/).

### Installation

Nach der Änderung der Quell-Datei müssen Sie die HTML-Datei `vortrag.html` erneut generieren. Hierzu können Sie wie folgt vorgehen:

1. Installieren Sie eine Python-Distribution. In diesem Beispiel verwenden wir Anaconda. Laden Sie also den Installer von [https://www.anaconda.com/download](https://www.anaconda.com/download) herunter und führen Sie den Installer aus.
2. Starten Sie vom Startmenü `Anaconda Prompt`. Erstellen Sie darüber eine neue virtuelle Umgebung, zum Beispiel mit dem Namen `quarto` wie im folgenden Beispiel, in der Sie später die notwendigen Pakete installieren können.

```bash
conda create --name quarto python
```

3. Installieren Sie Quarto. Laden Sie den Installer von [https://quarto.org/docs/get-started/](https://quarto.org/docs/get-started/) herunter und führen Sie den Installer aus. Installieren Sie NICHT das Pythonpaket `quarto-cli` via Python package manager pip!
4. Installieren Sie nun die notwendigen Python-Pakete, indem Sie vom Startmenü `Anaconda Prompt` starten und darin die folgenden Befehle ausführen:

```bash
conda activate quarto
pip install -r requirements.txt
```

Nun haben Sie eine fertige Umgebung, in der Sie die Quell-Datei `vortrag.qmd` editieren und veröffentlichen können. Statt der oben genannten manuellen Befehle können Sie auf Windows auch die Batch-Dateien `create_env.bat` verwenden. Achten Sie nur darauf davor den Schritt 3 auszuführen und die Batchdatei wie unter Schritt 2 beschrieben in einem Anaconda Prompt auszuführen.

### Veröffentlichen

Um die HTML-Datei `vortrag.html` zu veröffentlichen, führen Sie den folgenden Befehl in einem Anaconda Prompt aus (ersetzen Sie gegebenenfalls `quarto` durch den Namen Ihrer Anaconda-Umgebung):

```bash
conda activate quarto
quarto render vortrag.qmd
```

Alternativ können Sie auch die Batch-Datei `run_publish.bat` verwenden. Passen Sie in Zeile 6 gegebenenfalls das Installationsverzeichnis an. Diese Batch-Datei kann dann direkt per Doppelklick gestartet werden, ohne vorher einen Anaconda Prompt zu öffnen.

## Editieren der Simulatoren

Ergänzend zu dem Vortrag ist ein Simulator in zwei Varianten im Vortrag enthalten. Dieser zeigt beispielhafte Ergebnisse für Feder-Masse-Schwinger nach dem Lehrbuch Regelungstechnik 1 von Lunze. Die entsprechenden Dateien `simulation1.html` und `simulation2.html` sowie die zugehörigen JavaScript-Dateien `simulation1.js` und `simulation2.js` können mit einem beliebigen Texteditor wie [Visual Studio Code](https://code.visualstudio.com/) oder [Notepad++](https://notepad-plus-plus.org/) editiert werden.
