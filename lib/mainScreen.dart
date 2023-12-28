import 'dart:io';
import 'package:flutter/material.dart';
import 'package:image_picker/image_picker.dart';
import 'package:tflite/tflite.dart';
import 'package:percent_indicator/percent_indicator.dart';


class MyHomePage extends StatefulWidget {
  const MyHomePage({super.key, required this.title});

  final String title;

  @override
  State<MyHomePage> createState() => _MyHomePageState();
}

class _MyHomePageState extends State<MyHomePage> with TickerProviderStateMixin{
  File? _image;
  List? _outputs;
  bool loading = false;

  late AnimationController _animationController;
  late Animation<double> _animation;

  @override
  void initState() {
    super.initState();
    loading = true;

    _animationController = AnimationController(
      duration: Duration(seconds: 2),
      vsync: this,
    );

    _animation = Tween<double>(begin: 0.0, end: 1.0).animate(
      CurvedAnimation(
        parent: _animationController,
        curve: Curves.easeInOut,
      ),
    );

    _animationController.forward();


    loadModel().then((value) {
      setState(() {
        loading = false;
      });
    });
  }

  // Load the TFLite model
  Future<void> loadModel() async {
    await Tflite.loadModel(
        model: 'assets/brain-tumor-w.tflite',
        labels: 'assets/labels.txt',
        numThreads: 1);
  }

  //Run inference on an image
  runModelOnImage(File image) async {
    var output = await Tflite.runModelOnImage(
      path: image.path,
      imageMean: 0.0,
      imageStd: 255.0,
      numResults: 4,
      // Number of classification results
      threshold: 0.2,
      asynch: true,
    );
    setState(() {
      loading = false;
      _outputs = output;
    });
  }

  getImage() async {
    var pickedImage =
    await ImagePicker().pickImage(source: ImageSource.gallery);

    if (pickedImage == null) return null;

    setState(() {
      loading = true;
      _image = File(pickedImage.path);
    });
    runModelOnImage(_image!);
  }

  getImageFromCamera() async {
    var pickedImage = await ImagePicker().pickImage(source: ImageSource.camera);

    if (pickedImage == null) return null;

    setState(() {
      loading = true;
      _image = File(pickedImage.path);
    });
    runModelOnImage(_image!);
  }

  @override
  void dispose() {
    Tflite.close();
    super.dispose();
  }

  double extractFirstTwoDigits(double value) {
    // Convert the double to a string
    String stringValue = value.toString();

    // Find the index of the decimal point
    int decimalIndex = stringValue.indexOf('.');

    // Check if a decimal point exists and if there are at least two digits after it
    if (decimalIndex != -1 && decimalIndex + 3 <= stringValue.length) {
      // Extract the substring with the first two digits after the decimal point
      String extractedDigits = stringValue.substring(decimalIndex + 1, decimalIndex + 3);

      // Convert the extracted substring back to a double
      double result = double.parse('0.$extractedDigits');

      return result;
    }

    // Return the original value if there is no decimal point or not enough digits after it
    return value;
  }


  @override
  Widget build(BuildContext context) {
    return Scaffold(
      body: SingleChildScrollView(
        physics: BouncingScrollPhysics(),
        child: FadeTransition(
          opacity: _animation,
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.center,
            mainAxisAlignment: MainAxisAlignment.center,
            mainAxisSize: MainAxisSize.max,
            children: [
              Container(
                margin: EdgeInsets.only(top:20),
                padding: EdgeInsets.only(top:20),
                child: Image.asset(
                  'assets/logo.png',
                  height: 120,
                ),
              ), // Logo

              Container(
                decoration: BoxDecoration(
                    color: Colors.grey.shade200,
                    borderRadius: BorderRadius.all(Radius.circular(15))),
                clipBehavior: Clip.antiAliasWithSaveLayer,
                margin: EdgeInsets.only(right: 30, left: 30, bottom: 30),
                height: 300,
                width: double.infinity,
                child: _image != null
                    ? Image.file(
                  _image!,
                  // height: double.infinity,
                  // width: double.infinity,
                  alignment: Alignment.center,
                  fit: BoxFit.cover,
                )
                    : Container(
                    alignment: Alignment.center,
                    child: Text(
                      'Uploaded Image Will Appear Here.',
                      style: TextStyle(color: Colors.grey),
                    )),
              ), // Uploaded Photo

              Container(
                // padding: EdgeInsets.all(40),
                child: Row(
                  mainAxisAlignment: MainAxisAlignment.center,
                  crossAxisAlignment: CrossAxisAlignment.center,
                  children: [
                    Container(
                      margin: EdgeInsets.all(10),
                      child: FloatingActionButton(
                        onPressed: () {
                          getImageFromCamera();
                        },
                        child: Icon(
                          Icons.add_a_photo_rounded,
                          size: 30,
                        ),
                        backgroundColor: Colors.red.shade200,
                      ),
                    ),
                    Container(
                      margin: EdgeInsets.all(10),
                      child: FloatingActionButton(
                        onPressed: () {
                          getImage();
                        },
                        child: Icon(
                          Icons.add_photo_alternate_rounded,
                          size: 30,
                        ),
                        backgroundColor: Colors.red.shade200,
                      ),
                    ),
                  ],
                ),
              ), // Buttons
              Container(
                margin: EdgeInsets.all(20),
                padding: EdgeInsets.all(20),
                child: Column(
                  children: [
                    const Text(
                      'Results :',
                      style: TextStyle(
                          fontSize: 20,
                          fontWeight: FontWeight.bold,
                          color: Colors.grey),
                    ),
                    loading
                        ? Container(
                      margin: EdgeInsets.all(20),
                        padding: EdgeInsets.all(20),
                        child: CircularProgressIndicator(color: Colors.red))
                        : (_outputs != null && _outputs!.isNotEmpty)
                        ? Container(
                      padding: EdgeInsets.all(15),
                      child: Text(
                        '${_outputs![0]["label"]}',
                        style: TextStyle(
                          fontSize: 20,
                          color: Colors.red.shade900,
                          fontWeight: FontWeight.bold,
                        ),
                      ),
                    )
                        : const Text('No results yet.'),

                    (_outputs != null && _outputs!.isNotEmpty)
                        ? CircularPercentIndicator(
                      radius: 60.0,
                      circularStrokeCap: CircularStrokeCap.round,
                      animation: true,
                      percent: _outputs![0]['confidence'],
                      center: Text("${extractFirstTwoDigits(_outputs![0]['confidence']) * 100} %"),
                      progressColor: Colors.red,
                    )
                        : Text('')
                  ],
                ),
              ) // Result
            ],
          ),
        ),
      ),
    );
  }
}
