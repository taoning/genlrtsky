package main

import (
	"bytes"
	"encoding/binary"
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"text/template"
)

const (
	NSSAMP           = 20
	WVLSPAN          = 400 // Wavelength span in nm 780 - 380
	XRES             = 32
	YRES             = 32
	PROGRESSBARWIDTH = 50
	TOTALSTEPS       = YRES * XRES
)

type SimulationControl struct {
	OutName    string
	Version    bool
	NumThreads int
	Quiet      bool
}

type InputParams struct {
	LibPath             string
	SolFile             string
	GroundAlbedo        float64
	SolarZenith         string
	SolarAzimuth        string
	AerosolProfile      string
	AerosolOpticalDepth float64
	AodOverride         bool
	CosTheta            float64
	Phi                 float64
	CloudCover          float64
	CloudFile           string
	WavelengthGridFile  string
	Output              string
	Quiet               bool
}

type Datetime struct {
	Year   int
	Month  int
	Day    int
	Hour   int
	Minute int
}

type Location struct {
	Latitude         float64
	Longitude        float64
	StandardMeridian int
}

type Result struct {
	X, Y int
	Data [NSSAMP + 1]uint8 // Data to be written to hsrf
	Err  error
}

func panicError(e error) {
	if e != nil {
		panic(e)
	}
}

func handleError(e error) {
	if e != nil {
		fmt.Println(e)
		os.Exit(1)
	}
}

func parseMonth(str string) (int, error) {
	month, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if month < 1 || month > 12 {
		return 0, errors.New("Invalid month")
	}
	return month, nil
}

func parseDay(str string) (int, error) {
	day, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if day < 1 || day > 31 {
		return 0, errors.New("Invalid day")
	}
	return day, nil
}

func parseHour(str string) (int, error) {
	hours, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if hours < 0 || hours > 24 {
		return 0, errors.New("Invalid hours")
	}
	return hours, nil
}

func parseMinute(str string) (int, error) {
	minutes, err := strconv.Atoi(str)
	if err != nil {
		return 0, err
	}
	if minutes < 0 || minutes >= 60 {
		return 0, errors.New("Invalid minutes")
	}
	return minutes, nil
}

func parseAerosolProfile(input *InputParams) error {
	switch input.AerosolProfile {
	case "ca":
		input.AerosolProfile = "continental_average"
		return nil
	case "cc":
		input.AerosolProfile = "continental_clean"
		return nil
	case "cp":
		input.AerosolProfile = "continental_polluted"
		return nil
	case "mc":
		input.AerosolProfile = "maritime_clear"
		return nil
	case "mt":
		input.AerosolProfile = "maritime_tropical"
		return nil
	case "mp":
		input.AerosolProfile = "maritime_polluted"
		return nil
	case "a":
		input.AerosolProfile = "antarctic"
		return nil
	case "u":
		input.AerosolProfile = "urban"
		return nil
	case "d":
		input.AerosolProfile = "desert"
		return nil
	case "ds":
		input.AerosolProfile = "desert_spheroids"
		return nil
	default:
		return errors.New("Invalid aerosol profile, use one of these: ca, cc, cp, mc, mt, mp, a, u, d, ds")
	}
}

func printVersion() {
	fmt.Println("genlrtsky v0.0.2")

	cmd := exec.Command("uvspec", "-v")

	var stdoutBuf, stderrBuf bytes.Buffer
	cmd.Stdout = &stdoutBuf
	cmd.Stderr = &stderrBuf

	err := cmd.Run()
	panicError(err)
	fmt.Printf(stderrBuf.String())
	_, err = exec.LookPath("zenith")
	if err != nil {
		fmt.Println("zenith not found")
	}
	return
}

func getSolarAngle(input *InputParams, dt *Datetime, lc *Location) error {
	cmd := exec.Command("zenith", "-q", "-y", strconv.Itoa(dt.Year),
		"-a", strconv.FormatFloat(lc.Latitude, 'f', -1, 64),
		"-o", strconv.FormatFloat(lc.Longitude, 'f', -1, 64),
		"-s", strconv.Itoa(lc.StandardMeridian),
		strconv.Itoa(dt.Day), strconv.Itoa(dt.Month),
		strconv.Itoa(dt.Hour), strconv.Itoa(dt.Minute),
	)
	var stdoutBuf, stderrBuf bytes.Buffer
	cmd.Stdout = &stdoutBuf
	cmd.Stderr = &stderrBuf

	err := cmd.Run()
	fmt.Fprintf(os.Stderr, stderrBuf.String())
	solarAngles := strings.Fields(stdoutBuf.String())
	input.SolarZenith = solarAngles[1]
	input.SolarAzimuth = solarAngles[2]
	return err
}

func pix2loc(xres int, yres int, px int, py int) (float64, float64) {
	// Compute image location from pixel posistion
	// Assuming YMAJOR and YDECR
	x := px
	y := py
	y = yres - 1 - y
	loc0 := (float64(x) + 0.5) / float64(xres)
	loc1 := (float64(y) + 0.5) / float64(yres)
	return loc0, loc1
}

func viewray(x float64, y float64) (float64, float64) {
	vdir := [3]float64{0.0, 0.0, 1.0}
	hvec := [3]float64{-1.0, 0.0, 0.0}
	vvec := [3]float64{0.0, 1.0, 0.0}
	x -= 0.5
	y -= 0.5
	d := math.Sqrt(x*x + y*y)
	z := math.Cos(math.Pi * d)
	if d < 1e-6 {
		d = math.Pi
	} else {
		d = math.Sqrt(1-z*z) / d
	}
	x *= d
	y *= d
	direc := [3]float64{
		z*vdir[0] + x*hvec[0] + y*vvec[0],
		z*vdir[1] + x*hvec[1] + y*vvec[1],
		z*vdir[2] + x*hvec[2] + y*vvec[2],
	}
	phi := math.Atan2(direc[1], direc[0]) * 180 / math.Pi
	return -direc[2], phi
}

func scolor2scolr(scol [20]float64, ncs int) [21]uint8 {
	const COLXS uint8 = 128
	var p float64
	sclr := [21]uint8{}

	// Find the largest value in scol
	p = scol[0]
	for i := 1; i < ncs; i++ {
		if scol[i] > p {
			p = scol[i]
		}
	}
	if p > 1e-32 {
		frac, exp := math.Frexp(p)
		p = frac * 256.0 / p // Ensure this isn't zero by checking `p` and `frac`
		sclr[ncs] = uint8(exp) + COLXS
		for i := 0; i < ncs; i++ {
			if scol[i] > 0 {
				sclr[i] = uint8(scol[i] * p)
			}
		}
	} else {
		// Fill sclr with zeros if max value is not significant
		for i := range sclr {
			sclr[i] = 0
		}
	}
	return sclr
}

func worker(x, y int, input InputParams, tmpl *template.Template, wg *sync.WaitGroup, concurrencyLimits chan struct{}, results chan<- Result) {
	defer wg.Done()
	defer func() { <-concurrencyLimits }() // Release the "slot"

	var result Result

	result.X = x
	result.Y = y

	result.Data, result.Err = compute(x, y, input, tmpl)

	// Send the result to the results channel
	results <- result
}

func compute(x, y int, input InputParams, tmpl *template.Template) ([NSSAMP + 1]uint8, error) {

	loc0, loc1 := pix2loc(XRES, YRES, x, y)
	umu, phi := viewray(loc0, loc1)

	var templateBuffer bytes.Buffer
	input.CosTheta = umu
	input.Phi = phi
	err := tmpl.Execute(&templateBuffer, input)
	panicError(err)
	inputStr := templateBuffer.String()

	tdir, err := os.MkdirTemp("", "")
	panicError(err)
	defer os.RemoveAll(tdir)

	// fmt.Println(inputStr)
	cmd := exec.Command("uvspec")
	cmd.Stdin = bytes.NewBufferString(inputStr)
	cmd.Dir = tdir

	// Capture the output
	var stderr, out bytes.Buffer
	cmd.Stderr = &stderr
	cmd.Stdout = &out

	err = cmd.Run()
	panicError(err)
	res, err := os.ReadFile(filepath.Join(tdir, "mc.rad.spc"))
	panicError(err)
	rows := strings.Split(string(res), "\n")
	scolor := [20]float64{}
	for i, row := range rows {
		if row == "" {
			continue
		}
		fields := strings.Fields(row)
		value, err := strconv.ParseFloat(fields[4], 64)
		panicError(err)
		scolor[i] = value * WVLSPAN / 1000 // mW to W
	}
	sclr := scolor2scolr(scolor, NSSAMP)
	return sclr, err
}

func getSunDirection(input *InputParams) [3]float64 {
	sundir := [3]float64{}
	theta, err := strconv.ParseFloat(input.SolarZenith, 64)
	theta *= math.Pi / 180.0
	panicError(err)
	phi, err := strconv.ParseFloat(input.SolarAzimuth, 64)
	phi = (270 - phi) * math.Pi / 180.0
	panicError(err)
	sundir[0] = math.Sin(theta) * math.Cos(phi)
	sundir[1] = math.Sin(theta) * math.Sin(phi)
	sundir[2] = math.Cos(theta)
	return sundir
}

func parsePositionalArgs(datetime *Datetime, args []string) error {
	// Check and get positional arguments before any optional flags.
	positionalArgs := []string{}
	for _, arg := range args { // Skipping the program name.
		if strings.HasPrefix(arg, "-") {
			break // Stop on the first flag.
		}
		positionalArgs = append(positionalArgs, arg)
	}

	if len(positionalArgs) != 4 {
		flag.Usage()
		return errors.New("Invalid number of positional arguments")
	}
	month, err := parseMonth(positionalArgs[0])
	handleError(err)
	day, err := parseDay(positionalArgs[1])
	handleError(err)
	hour, err := parseHour(positionalArgs[2])
	handleError(err)
	minute, err := parseMinute(positionalArgs[3])
	handleError(err)
	datetime.Month = month
	datetime.Day = day
	datetime.Hour = hour
	datetime.Minute = minute
	return err
}

func setupFlags(simctrl *SimulationControl, input *InputParams, dateTime *Datetime, location *Location) {
	flag.IntVar(&dateTime.Year, "y", 2003, "Year")
	flag.Float64Var(&location.Latitude, "a", 37.0, "Latitude, north positive")
	flag.Float64Var(&location.Longitude, "o", 122.0, "Longitude, west positive")
	flag.IntVar(&location.StandardMeridian, "m", 120, "Standard meridian, west positive")
	flag.Float64Var(&input.GroundAlbedo, "g", 0.2, "Albedo")
	flag.StringVar(&input.AerosolProfile, "s", "ca", "Standard aerosol profile name")
	flag.Float64Var(&input.AerosolOpticalDepth, "d", 0.0, "Aerosol optical depth")
	flag.StringVar(&simctrl.OutName, "p", "out", "Output file prefix")
	flag.BoolVar(&simctrl.Version, "version", false, "Print version")
	flag.Float64Var(&input.CloudCover, "c", 0.0, "Cloud cover 0-1")
	flag.BoolVar(&simctrl.Quiet, "quiet", false, "Quiet mode")
	flag.IntVar(&simctrl.NumThreads, "n", 1, "Number of threads")

	// Customize the usage function
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: genlrtsky month day hour minute [options]\n")
		flag.PrintDefaults()
	}
}

func getLibPathOrExit() string {
	libPath := os.Getenv("LIBRADTRAN_DATA_FILES")
	if libPath == "" {
		fmt.Fprintln(os.Stderr, "LIBRADTRAN_DATA_FILES not set")
		os.Exit(1)
	}
	return libPath
}

func displayProgressBar(completedSteps int32, totalSteps int, progressBarWidth int) {
	progress := float64(completedSteps) / float64(totalSteps)
	bar := strings.Repeat("=", int(progress*float64(progressBarWidth)))
	fmt.Printf("\r[%-50s] %d%%", bar, int(progress*100))
}

func writeHsr(hsrf *os.File, hsr *bytes.Buffer) {
	hsrf.WriteString("#?RADIANCE\n")
	fmt.Fprintf(hsrf, "NCOMP=%d\n", NSSAMP)
	hsrf.WriteString("WAVELENGTH_SPLITS= 380 480 588 780\n")
	hsrf.WriteString("FORMAT=Radiance_spectra\n\n")
	fmt.Fprintf(hsrf, "-Y       %d +X       %d\n", YRES, XRES)
	hsrf.Write(hsr.Bytes())
}

func main() {

	completedSteps := int32(0)

	inputTemplate := `data_files_path {{.LibPath}} 
source solar {{.SolFile}}
rte_solver mystic
albedo {{.GroundAlbedo}}
sza {{.SolarZenith}}
phi0 {{.SolarAzimuth}}
aerosol_default
aerosol_species_file {{.AerosolProfile}}
{{if .AodOverride}}aerosol_modify tau set {{.AerosolOpticalDepth}}{{end}}
{{if .CloudCover}}
wc_file 1D {{.CloudFile}}
cloudcover wc {{.CloudCover}}
interpret_as_level wc
{{end}}
umu {{.CosTheta}}
phi {{.Phi}}
wavelength_grid_file {{.WavelengthGridFile}}
mol_abs_param crs
mc_spherical 1D
mc_vroom on
{{if .Quiet}}quiet{{end}}`

	tmpl, err := template.New("input").Parse(inputTemplate)
	panicError(err)

	input := InputParams{}

	// Extraterrestrial solar irradiance in mW/m^2/nm
	EXTSOL := [NSSAMP]float64{1192.74, 1772.2, 1394.52, 2069.67,
		2074.3, 1964.5, 1971.29, 1942.59, 1870.53,
		1852.63, 1772.84, 1737.29, 1662.2, 1586.51,
		1537.19, 1475.65, 1400.21, 1320.18, 1273.23,
		1220.45,
	}

	input.LibPath = getLibPathOrExit()
	dateTime := Datetime{}
	location := Location{}
	simctrl := SimulationControl{}

	setupFlags(&simctrl, &input, &dateTime, &location)

	if len(os.Args) < 2 {
		flag.Usage()
		os.Exit(0)
	}

	if os.Args[1] == "-version" {
		printVersion()
		os.Exit(0)
	}

	if len(os.Args) < 5 {
		flag.Usage()
		os.Exit(0)
	}

	positionalArgs := os.Args[1:5]
	os.Args = os.Args[4:]
	flag.Parse()
	if err = parsePositionalArgs(&dateTime, positionalArgs); err != nil {
		handleError(err)
	}

	err = parseAerosolProfile(&input)
	handleError(err)
	input.AodOverride = false
	if input.AerosolOpticalDepth > 0 {
		input.AodOverride = true
	}

	if simctrl.Quiet {
		input.Quiet = true
	}

	outRad := simctrl.OutName + ".rad"
	outHsr := simctrl.OutName + ".hsr"

	hsrf, err := os.Create(outHsr)
	handleError(err)
	defer hsrf.Close()

	radf, err := os.Create(outRad)
	handleError(err)
	defer radf.Close()

	// Get the solar angle
	err = getSolarAngle(&input, &dateTime, &location)
	panicError(err)

	// Create a temp file to store solar and wavelength data
	wavelengthGridFile, err := os.CreateTemp("", "lambda")
	panicError(err)
	defer os.Remove(wavelengthGridFile.Name())
	defer wavelengthGridFile.Close()
	input.WavelengthGridFile = wavelengthGridFile.Name()

	solFile, err := os.CreateTemp("", "solar")
	panicError(err)
	defer os.Remove(solFile.Name())
	defer solFile.Close()
	input.SolFile = solFile.Name()

	for index, value := range EXTSOL {
		wavelength := 390 + index*NSSAMP
		wavelengthStr := strconv.Itoa(wavelength)
		wavelengthGridFile.WriteString(wavelengthStr + "\n")
		solFile.WriteString(wavelengthStr + " " + strconv.FormatFloat(value, 'f', -1, 64) + "\n")
	}

	if input.CloudCover > 0 {
		cloudFile, err := os.CreateTemp("", "cloud")
		panicError(err)
		defer os.Remove(cloudFile.Name())
		defer cloudFile.Close()
		cloudFile.WriteString(`
		5.0	0	0
		4.0	0.	0.0
		3.0	0.25	10.0
		2.0	0.	0.0
		1.0	0	0.0
		0.0	0	0.0`)
		input.CloudFile = cloudFile.Name()
	}
	results := make(chan Result, YRES*XRES)
	wg := &sync.WaitGroup{}

	// Limit the number of concurrent goroutines to the number of processor cores
	numCores := min(simctrl.NumThreads, runtime.NumCPU())
	concurrencyLimit := make(chan struct{}, numCores) // Semaphore-like channel

	for y := 0; y < YRES; y++ {
		for x := 0; x < XRES; x++ {
			concurrencyLimit <- struct{}{} // Acquire a "slot"
			wg.Add(1)
			go worker(x, y, input, tmpl, wg, concurrencyLimit, results)
			atomic.AddInt32(&completedSteps, 1)
			if !simctrl.Quiet {
				displayProgressBar(completedSteps, TOTALSTEPS, PROGRESSBARWIDTH)
			}
		}
	}

	// Ensure all workers are finished before closing the results channel
	go func() {
		wg.Wait()
		close(results)
	}()

	orderedResults := make([][]Result, YRES)
	for y := range orderedResults {
		orderedResults[y] = make([]Result, XRES)
	}

	for res := range results {
		if res.Err != nil {
			fmt.Printf("Error processing (%d, %d): %v\n", res.X, res.Y, res.Err)
			continue
		}
		orderedResults[res.Y][res.X] = res
	}

	// Write the results in `hsrf` in order
	var hsr bytes.Buffer
	for y := range orderedResults {
		for x := range orderedResults[y] {
			if res := orderedResults[y][x]; res.Err == nil {
				for _, value := range res.Data {
					err := binary.Write(&hsr, binary.LittleEndian, value)
					panicError(err)
				}
			}
		}
	}

	var templateBuffer bytes.Buffer
	input.CosTheta = -1
	input.Phi = 0
	err = tmpl.Execute(&templateBuffer, input)
	panicError(err)
	inputStr := templateBuffer.String()
	inputStr += "\n" + "output_user edir"
	inputStr += "\n" + "rte_solver disort"

	tdir, err := os.MkdirTemp("", "")
	panicError(err)
	defer os.RemoveAll(tdir)

	// fmt.Println(inputStr)
	cmd := exec.Command("uvspec")
	cmd.Stdin = bytes.NewBufferString(inputStr)
	cmd.Dir = tdir

	// Capture the output
	var stderr, out bytes.Buffer
	cmd.Stderr = &stderr
	cmd.Stdout = &out

	err = cmd.Run()
	panicError(err)
	sun_radiance := [NSSAMP]float64{}
	for idx, value := range strings.Split(strings.Trim(out.String(), " \n"), "\n") {
		sun_radiance[idx], err = strconv.ParseFloat(strings.TrimSpace(value), 64)
		panicError(err)
	}

	sundir := getSunDirection(&input)
	writeHsr(hsrf, &hsr)

	fmt.Fprintf(radf, "# Latitude: %f, Longitude: %f, TimeZone: %d\n", location.Latitude, location.Longitude, location.StandardMeridian)
	fmt.Fprintf(radf, "# Date time: %d-%d-%d %d:%02d\n", dateTime.Year, dateTime.Month, dateTime.Day, dateTime.Hour, dateTime.Minute)
	fmt.Fprintf(radf, "# Solar zenith and azimuth: %s, %s\n", input.SolarZenith, input.SolarAzimuth)
	fmt.Fprintf(radf, "# Ground albedo: %f\n", input.GroundAlbedo)
	fmt.Fprintf(radf, "# Aerosol profile: %s\n", input.AerosolProfile)
	if input.AodOverride {
		fmt.Fprintf(radf, "# Overriding aerosol optical depth: %f\n", input.AerosolOpticalDepth)
	}
	if input.CloudCover > 0 {
		fmt.Fprintf(radf, "# Cloud cover: %f\n", input.CloudCover)
	}
	radf.WriteString("void spectrum sun\n0\n0\n22 380 780 ")
	for _, value := range sun_radiance {
		fmt.Fprintf(radf, "%.1f ", value*WVLSPAN/6.7967e-5/1000/sundir[2])
	}
	radf.WriteString("\n\nsun light solar\n0\n0\n3 1 1 1\n\n")
	fmt.Fprintf(radf, "solar source sun\n0\n0\n4 %f %f %f 0.533\n\n", sundir[0], sundir[1], sundir[2])
	fmt.Fprintf(radf, "void specpict skyfunc 9 noop %s fisheye.cal fish_u fish_v -rx 90 -rz 90 0 0\n", outHsr)
	radf.WriteString("skyfunc glow skyglow 0 0 4 1 1 1 0\n")
	radf.WriteString("skyglow source sky 0 0 4 0 0 1 180\n")
}
